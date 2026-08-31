!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Restart file for k point calculations
!> \par History
!> \author JGH (30.09.2015)
! **************************************************************************************************
MODULE kpoint_io

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_release,&
                                              cp_fm_type,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              default_string_length
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: nso
   USE particle_types,                  ONLY: particle_type
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_io,                        ONLY: wfn_restart_file_name
   USE qs_ot_types,                     ONLY: qs_ot_check_channel_context,&
                                              qs_ot_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              qs_scf_ot_restart_clear
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_io'

   INTEGER, PARAMETER                   :: version = 1
   INTEGER, PARAMETER                   :: ot_restart_version = 1
   CHARACTER(LEN=*), PARAMETER          :: ot_restart_marker = "CP2K_KPOINT_OT_RESTART"

   PUBLIC :: read_kpoints_restart, &
             write_kpoints_restart
   PUBLIC :: get_cell, write_kpoints_file_header

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param denmat ...
!> \param kpoints ...
!> \param scf_env ...
!> \param dft_section ...
!> \param particle_set ...
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE write_kpoints_restart(denmat, kpoints, scf_env, dft_section, &
                                    particle_set, qs_kind_set)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: denmat
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_kpoints_restart'
      CHARACTER(LEN=30), DIMENSION(2), PARAMETER :: &
         keys = ["SCF%PRINT%RESTART_HISTORY", "SCF%PRINT%RESTART        "]

      INTEGER                                            :: handle, ic, ikey, ires, ispin, nimages, &
                                                            nspin
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(cp_fm_type), POINTER                          :: fmat
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           dft_section, keys(1)), cp_p_file) .OR. &
          BTEST(cp_print_key_should_output(logger%iter_info, &
                                           dft_section, keys(2)), cp_p_file)) THEN

         DO ikey = 1, SIZE(keys)

            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 dft_section, keys(ikey)), cp_p_file)) THEN

               ires = cp_print_key_unit_nr(logger, dft_section, keys(ikey), &
                                           extension=".kp", file_status="REPLACE", file_action="WRITE", &
                                           do_backup=.TRUE., file_form="UNFORMATTED")

               CALL write_kpoints_file_header(qs_kind_set, particle_set, ires)

               nspin = SIZE(denmat, 1)
               nimages = SIZE(denmat, 2)
               NULLIFY (cell_to_index)
               IF (nimages > 1) THEN
                  CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
               END IF
               CPASSERT(ASSOCIATED(scf_env%scf_work1))
               fmat => scf_env%scf_work1(1)

               DO ispin = 1, nspin
                  IF (ires > 0) WRITE (ires) ispin, nspin, nimages
                  DO ic = 1, nimages
                     IF (nimages > 1) THEN
                        cell = get_cell(ic, cell_to_index)
                     ELSE
                        cell = 0
                     END IF
                     IF (ires > 0) WRITE (ires) ic, cell
                     CALL copy_dbcsr_to_fm(denmat(ispin, ic)%matrix, fmat)
                     CALL cp_fm_write_unformatted(fmat, ires)
                  END DO
               END DO
               IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
                  CALL write_kpoint_ot_restart_payload(scf_env%qs_ot_env, ires)
               END IF

               CALL cp_print_key_finished_output(ires, logger, dft_section, TRIM(keys(ikey)))

            END IF

         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE write_kpoints_restart

! **************************************************************************************************
!> \brief ...
!> \param ic ...
!> \param cell_to_index ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_cell(ic, cell_to_index) RESULT(cell)
      INTEGER, INTENT(IN)                                :: ic
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, DIMENSION(3)                              :: cell

      INTEGER                                            :: i1, i2, i3

      cell = 0
      ALLCELL: DO i3 = LBOUND(cell_to_index, 3), UBOUND(cell_to_index, 3)
         DO i2 = LBOUND(cell_to_index, 2), UBOUND(cell_to_index, 2)
            DO i1 = LBOUND(cell_to_index, 1), UBOUND(cell_to_index, 1)
               IF (cell_to_index(i1, i2, i3) == ic) THEN
                  cell(1) = i1
                  cell(2) = i2
                  cell(3) = i3
                  EXIT ALLCELL
               END IF
            END DO
         END DO
      END DO ALLCELL

   END FUNCTION get_cell

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param ires ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE write_kpoints_file_header(qs_kind_set, particle_set, ires, basis_type)

      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER                                            :: ires
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_kpoints_file_header'

      CHARACTER(LEN=default_string_length)               :: basis_type_local
      INTEGER                                            :: handle, iatom, ikind, iset, ishell, &
                                                            lmax, lshell, nao, natom, nset, &
                                                            nset_max, nsgf, nshell_max
      INTEGER, DIMENSION(:), POINTER                     :: nset_info, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l, nshell_info
      INTEGER, DIMENSION(:, :, :), POINTER               :: nso_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_parameter

      CALL timeset(routineN, handle)

      IF (ires > 0) THEN

         IF (PRESENT(basis_type)) THEN
            basis_type_local = basis_type
         ELSE
            basis_type_local = "ORB"
         END IF

         ! create some info about the basis set first
         natom = SIZE(particle_set, 1)
         nset_max = 0
         nshell_max = 0

         DO iatom = 1, natom
            NULLIFY (orb_basis_set, dftb_parameter)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_type=basis_type_local, &
                             basis_set=orb_basis_set, dftb_parameter=dftb_parameter)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      nshell=nshell, &
                                      l=l)
               nset_max = MAX(nset_max, nset)
               DO iset = 1, nset
                  nshell_max = MAX(nshell_max, nshell(iset))
               END DO
            ELSE IF (ASSOCIATED(dftb_parameter)) THEN
               CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
               nset_max = MAX(nset_max, 1)
               nshell_max = MAX(nshell_max, lmax + 1)
            ELSE
               CPABORT("Unknown basis type.")
            END IF
         END DO

         ALLOCATE (nso_info(nshell_max, nset_max, natom))
         nso_info(:, :, :) = 0

         ALLOCATE (nshell_info(nset_max, natom))
         nshell_info(:, :) = 0

         ALLOCATE (nset_info(natom))
         nset_info(:) = 0

         nao = 0
         DO iatom = 1, natom
            NULLIFY (orb_basis_set, dftb_parameter)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_type=basis_type_local, &
                             basis_set=orb_basis_set, dftb_parameter=dftb_parameter)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, nsgf=nsgf, nset=nset, nshell=nshell, l=l)
               nao = nao + nsgf
               nset_info(iatom) = nset
               DO iset = 1, nset
                  nshell_info(iset, iatom) = nshell(iset)
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     nso_info(ishell, iset, iatom) = nso(lshell)
                  END DO
               END DO
            ELSE IF (ASSOCIATED(dftb_parameter)) THEN
               CALL get_dftb_atom_param(dftb_parameter, lmax=lmax)
               nset_info(iatom) = 1
               nshell_info(1, iatom) = lmax + 1
               DO ishell = 1, lmax + 1
                  lshell = ishell - 1
                  nso_info(ishell, 1, iatom) = nso(lshell)
               END DO
               nao = nao + (lmax + 1)**2
            ELSE
               CPABORT("Unknown basis set type.")
            END IF
         END DO

         WRITE (ires) version
         WRITE (ires) natom, nao, nset_max, nshell_max
         WRITE (ires) nset_info
         WRITE (ires) nshell_info
         WRITE (ires) nso_info

         DEALLOCATE (nset_info, nshell_info, nso_info)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_kpoints_file_header

! **************************************************************************************************
!> \brief ...
!> \param denmat ...
!> \param kpoints ...
!> \param fmwork ...
!> \param natom ...
!> \param para_env ...
!> \param id_nr ...
!> \param dft_section ...
!> \param natom_mismatch ...
! **************************************************************************************************
   SUBROUTINE read_kpoints_restart(denmat, kpoints, fmwork, natom, &
                                   para_env, id_nr, dft_section, natom_mismatch, scf_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: denmat
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: fmwork
      INTEGER, INTENT(IN)                                :: natom
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: id_nr
      TYPE(section_vals_type), POINTER                   :: dft_section
      LOGICAL, INTENT(OUT)                               :: natom_mismatch
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_kpoints_restart'

      CHARACTER(LEN=default_path_length)                 :: file_name
      INTEGER                                            :: handle, restart_unit
      LOGICAL                                            :: exist
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      restart_unit = -1
      IF (PRESENT(scf_env)) CALL qs_scf_ot_restart_clear(scf_env)

      IF (para_env%is_source()) THEN

         CALL wfn_restart_file_name(file_name, exist, dft_section, logger, kp=.TRUE.)
         IF (id_nr /= 0) THEN
            ! Is it one of the backup files?
            file_name = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(id_nr))
         END IF

         CALL open_file(file_name=file_name, &
                        file_action="READ", &
                        file_form="UNFORMATTED", &
                        file_status="OLD", &
                        unit_number=restart_unit)

      END IF

      IF (PRESENT(scf_env)) THEN
         CALL read_kpoints_restart_low(denmat, kpoints, fmwork(1), para_env, &
                                       natom, restart_unit, natom_mismatch, scf_env=scf_env)
      ELSE
         CALL read_kpoints_restart_low(denmat, kpoints, fmwork(1), para_env, &
                                       natom, restart_unit, natom_mismatch)
      END IF

      ! only the io_node returns natom_mismatch, must broadcast it
      CALL para_env%bcast(natom_mismatch)

      ! Close restart file
      IF (para_env%is_source()) CALL close_file(unit_number=restart_unit)

      CALL timestop(handle)

   END SUBROUTINE read_kpoints_restart

! **************************************************************************************************
!> \brief Reading the mos from apreviously defined restart file
!> \param denmat ...
!> \param kpoints ...
!> \param fmat ...
!> \param para_env ...
!> \param natom ...
!> \param rst_unit ...
!> \param natom_mismatch ...
!> \par History
!>      12.2007 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE read_kpoints_restart_low(denmat, kpoints, fmat, para_env, &
                                       natom, rst_unit, natom_mismatch, scf_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: denmat
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fmat
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: natom, rst_unit
      LOGICAL, INTENT(OUT)                               :: natom_mismatch
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env

      INTEGER                                            :: ic, image, ispin, nao, nao_read, &
                                                            natom_read, ni, nimages, nset_max, &
                                                            nshell_max, nspin, nspin_read, &
                                                            version_read
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: natom_match

      CPASSERT(ASSOCIATED(denmat))
      CALL dbcsr_get_info(denmat(1, 1)%matrix, nfullrows_total=nao)

      nspin = SIZE(denmat, 1)
      nimages = SIZE(denmat, 2)
      NULLIFY (cell_to_index)
      IF (nimages > 1) THEN
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      IF (para_env%is_source()) THEN
         READ (rst_unit) version_read
         IF (version_read /= version) THEN
            CALL cp_abort(__LOCATION__, &
                          " READ RESTART : Version of restart file not supported")
         END IF
         READ (rst_unit) natom_read, nao_read, nset_max, nshell_max
         natom_match = (natom_read == natom)
         IF (natom_match) THEN
            IF (nao_read /= nao) THEN
               CALL cp_abort(__LOCATION__, &
                             " READ RESTART : Different number of basis functions")
            END IF
            READ (rst_unit)
            READ (rst_unit)
            READ (rst_unit)
         END IF
      END IF

      ! make natom_match and natom_mismatch uniform across all nodes
      CALL para_env%bcast(natom_match)
      natom_mismatch = .NOT. natom_match
      ! handle natom_match false
      IF (natom_mismatch) THEN
         CALL cp_warn(__LOCATION__, &
                      " READ RESTART : WARNING : DIFFERENT natom, returning ")
      ELSE

         DO
            ispin = 0
            IF (para_env%is_source()) THEN
               READ (rst_unit) ispin, nspin_read, ni
            END IF
            CALL para_env%bcast(ispin)
            CALL para_env%bcast(nspin_read)
            CALL para_env%bcast(ni)
            IF (nspin_read /= nspin) THEN
               CALL cp_abort(__LOCATION__, &
                             " READ RESTART : Different spin polarisation not supported")
            END IF
            DO ic = 1, ni
               IF (para_env%is_source()) THEN
                  READ (rst_unit) image, cell
               END IF
               CALL para_env%bcast(image)
               CALL para_env%bcast(cell)
               !
               CALL cp_fm_read_unformatted(fmat, rst_unit)
               !
               IF (nimages > 1) THEN
                  image = get_index(cell, cell_to_index)
               ELSE IF (SUM(ABS(cell(:))) == 0) THEN
                  image = 1
               ELSE
                  image = 0
               END IF
               IF (image >= 1 .AND. image <= nimages) THEN
                  CALL copy_fm_to_dbcsr(fmat, denmat(ispin, image)%matrix)
               END IF
            END DO
            IF (ispin == nspin_read) EXIT
         END DO
         IF (PRESENT(scf_env)) THEN
            CALL read_kpoint_ot_restart_payload(rst_unit, para_env, fmat, scf_env)
         END IF

      END IF

   END SUBROUTINE read_kpoints_restart_low

! **************************************************************************************************
!> \brief write optional k-point OT restart state after the density payload
!> \param qs_ot_env OT channels
!> \param ires restart unit
! **************************************************************************************************
   SUBROUTINE write_kpoint_ot_restart_payload(qs_ot_env, ires)
      TYPE(qs_ot_type), DIMENSION(:), INTENT(IN)         :: qs_ot_env
      INTEGER, INTENT(IN)                                :: ires

      INTEGER                                            :: ichannel, nchannel, nener, nmo, nmo_ortho
      INTEGER, DIMENSION(12)                             :: int_payload
      LOGICAL                                            :: has_ener_x, has_rot_mat_x, &
                                                            has_rot_mat_x_im
      LOGICAL, DIMENSION(9)                              :: log_payload
      REAL(KIND=dp), DIMENSION(8)                        :: real_payload

      nchannel = SIZE(qs_ot_env)
      IF (ires > 0) THEN
         WRITE (ires) ot_restart_marker
         WRITE (ires) ot_restart_version, nchannel
      END IF
      DO ichannel = 1, nchannel
         nmo = 0
         nmo_ortho = 0
         IF (ASSOCIATED(qs_ot_env(ichannel)%matrix_c0)) THEN
            CALL dbcsr_get_info(qs_ot_env(ichannel)%matrix_c0, nfullcols_total=nmo)
         END IF
         IF (ASSOCIATED(qs_ot_env(ichannel)%matrix_sc0)) THEN
            CALL dbcsr_get_info(qs_ot_env(ichannel)%matrix_sc0, nfullcols_total=nmo_ortho)
         END IF
         has_rot_mat_x = qs_ot_env(ichannel)%state_allocated .AND. &
                         ASSOCIATED(qs_ot_env(ichannel)%rot_mat_x)
         has_rot_mat_x_im = qs_ot_env(ichannel)%state_allocated .AND. &
                            ASSOCIATED(qs_ot_env(ichannel)%rot_mat_x_im)
         has_ener_x = qs_ot_env(ichannel)%state_allocated .AND. ASSOCIATED(qs_ot_env(ichannel)%ener_x)
         nener = 0
         IF (has_ener_x) nener = SIZE(qs_ot_env(ichannel)%ener_x)

         int_payload = [ichannel, qs_ot_env(ichannel)%spin_index, &
                        qs_ot_env(ichannel)%kpoint_index, &
                        qs_ot_env(ichannel)%local_kpoint_index, nmo, nmo_ortho, &
                        qs_ot_env(ichannel)%line_search_left, &
                        qs_ot_env(ichannel)%line_search_right, &
                        qs_ot_env(ichannel)%line_search_mid, &
                        qs_ot_env(ichannel)%line_search_count, &
                        qs_ot_env(ichannel)%diis_iter, nener]
         log_payload = [qs_ot_env(ichannel)%state_allocated, &
                        qs_ot_env(ichannel)%has_complex_kpoint_state, &
                        qs_ot_env(ichannel)%settings%do_rotation, &
                        qs_ot_env(ichannel)%settings%do_ener, &
                        has_rot_mat_x, has_rot_mat_x_im, has_ener_x, &
                        qs_ot_env(ichannel)%line_search_might_be_done, &
                        qs_ot_env(ichannel)%energy_only]
         real_payload = [qs_ot_env(ichannel)%kpoint_weight, &
                         qs_ot_env(ichannel)%delta, qs_ot_env(ichannel)%gnorm, &
                         qs_ot_env(ichannel)%gnorm_old, qs_ot_env(ichannel)%etotal, &
                         qs_ot_env(ichannel)%gradient, qs_ot_env(ichannel)%ds_min, &
                         qs_ot_env(ichannel)%broyden_adaptive_sigma]

         IF (ires > 0) THEN
            WRITE (ires) int_payload
            WRITE (ires) log_payload
            WRITE (ires) real_payload
            WRITE (ires) qs_ot_env(ichannel)%ot_pos
            WRITE (ires) qs_ot_env(ichannel)%ot_energy
            WRITE (ires) qs_ot_env(ichannel)%ot_grad
         END IF
         IF (has_rot_mat_x) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_x, &
                                                            qs_ot_env(ichannel), ires)
         IF (has_rot_mat_x_im) CALL write_restart_dbcsr_matrix(qs_ot_env(ichannel)%rot_mat_x_im, &
                                                               qs_ot_env(ichannel), ires)
         IF (has_ener_x .AND. ires > 0) WRITE (ires) qs_ot_env(ichannel)%ener_x
      END DO

   END SUBROUTINE write_kpoint_ot_restart_payload

! **************************************************************************************************
!> \brief write a DBCSR matrix through the restart FM layer
!> \param matrix matrix to write
!> \param qs_ot_env channel owning the matrix distribution
!> \param ires restart unit
! **************************************************************************************************
   SUBROUTINE write_restart_dbcsr_matrix(matrix, qs_ot_env, ires)
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      INTEGER, INTENT(IN)                                :: ires

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fmat

      NULLIFY (fm_struct)
      CALL dbcsr_get_info(matrix, nfullrows_total=nrow, nfullcols_total=ncol)
      CALL cp_fm_struct_create(fm_struct, para_env=qs_ot_env%para_env, &
                               context=qs_ot_env%blacs_env, &
                               nrow_global=nrow, ncol_global=ncol)
      CALL cp_fm_create(fmat, fm_struct)
      CALL copy_dbcsr_to_fm(matrix, fmat)
      CALL cp_fm_write_unformatted(fmat, ires)
      CALL cp_fm_release(fmat)
      CALL cp_fm_struct_release(fm_struct)

   END SUBROUTINE write_restart_dbcsr_matrix

! **************************************************************************************************
!> \brief read optional k-point OT restart state after the density payload
!> \param rst_unit restart unit
!> \param para_env parallel environment
!> \param fmat restart FM template
!> \param scf_env SCF environment receiving the cache
! **************************************************************************************************
   SUBROUTINE read_kpoint_ot_restart_payload(rst_unit, para_env, fmat, scf_env)
      INTEGER, INTENT(IN)                                :: rst_unit
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fmat
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(LEN=LEN(ot_restart_marker))              :: marker
      INTEGER                                            :: ichannel, ios, nchannel, nener, &
                                                            version_read
      INTEGER, DIMENSION(12)                             :: int_payload
      LOGICAL                                            :: payload_found
      LOGICAL, DIMENSION(9)                              :: log_payload
      REAL(KIND=dp), DIMENSION(8)                        :: real_payload

      payload_found = .FALSE.
      marker = ""
      IF (para_env%is_source()) THEN
         READ (rst_unit, IOSTAT=ios) marker
         payload_found = (ios == 0 .AND. TRIM(marker) == ot_restart_marker)
      END IF
      CALL para_env%bcast(payload_found)
      IF (.NOT. payload_found) RETURN

      version_read = 0
      nchannel = 0
      IF (para_env%is_source()) READ (rst_unit) version_read, nchannel
      CALL para_env%bcast(version_read)
      CALL para_env%bcast(nchannel)
      IF (version_read /= ot_restart_version) THEN
         CALL cp_warn(__LOCATION__, "Ignoring unsupported k-point OT restart payload version")
         RETURN
      END IF

      CALL qs_scf_ot_restart_clear(scf_env)
      scf_env%qs_ot_restart%valid = .TRUE.
      scf_env%qs_ot_restart%format_version = version_read
      scf_env%qs_ot_restart%nchannel = nchannel
      ALLOCATE (scf_env%qs_ot_restart%channel(nchannel))
      DO ichannel = 1, nchannel
         IF (para_env%is_source()) THEN
            READ (rst_unit) int_payload
            READ (rst_unit) log_payload
            READ (rst_unit) real_payload
            READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ot_pos
            READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ot_energy
            READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ot_grad
         END IF
         CALL para_env%bcast(int_payload)
         CALL para_env%bcast(log_payload)
         CALL para_env%bcast(real_payload)
         CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ot_pos)
         CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ot_energy)
         CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ot_grad)

         scf_env%qs_ot_restart%channel(ichannel)%valid = log_payload(1)
         scf_env%qs_ot_restart%channel(ichannel)%spin_index = int_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%kpoint_index = int_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%local_kpoint_index = int_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%nmo = int_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%nmo_ortho = int_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_left = int_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_right = int_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_mid = int_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_count = int_payload(10)
         scf_env%qs_ot_restart%channel(ichannel)%diis_iter = int_payload(11)
         nener = int_payload(12)
         scf_env%qs_ot_restart%channel(ichannel)%has_complex_kpoint_state = log_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%do_rotation = log_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%do_ener = log_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x = log_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x_im = log_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%has_ener_x = log_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%line_search_might_be_done = log_payload(8)
         scf_env%qs_ot_restart%channel(ichannel)%energy_only = log_payload(9)
         scf_env%qs_ot_restart%channel(ichannel)%kpoint_weight = real_payload(1)
         scf_env%qs_ot_restart%channel(ichannel)%delta = real_payload(2)
         scf_env%qs_ot_restart%channel(ichannel)%gnorm = real_payload(3)
         scf_env%qs_ot_restart%channel(ichannel)%gnorm_old = real_payload(4)
         scf_env%qs_ot_restart%channel(ichannel)%etotal = real_payload(5)
         scf_env%qs_ot_restart%channel(ichannel)%gradient = real_payload(6)
         scf_env%qs_ot_restart%channel(ichannel)%ds_min = real_payload(7)
         scf_env%qs_ot_restart%channel(ichannel)%broyden_adaptive_sigma = real_payload(8)

         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x) THEN
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%nmo, &
                                        scf_env%qs_ot_restart%channel(ichannel)%nmo, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x, rst_unit)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_rot_mat_x_im) THEN
            CALL read_restart_fm_matrix(fmat, scf_env%qs_ot_restart%channel(ichannel)%nmo, &
                                        scf_env%qs_ot_restart%channel(ichannel)%nmo, &
                                        scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_im, rst_unit)
         END IF
         IF (scf_env%qs_ot_restart%channel(ichannel)%has_ener_x) THEN
            ALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_x(nener))
            IF (para_env%is_source()) READ (rst_unit) scf_env%qs_ot_restart%channel(ichannel)%ener_x
            CALL para_env%bcast(scf_env%qs_ot_restart%channel(ichannel)%ener_x)
         END IF
      END DO

   END SUBROUTINE read_kpoint_ot_restart_payload

! **************************************************************************************************
!> \brief read a restart FM matrix using an existing FM distribution as template
!> \param fmat_template FM template
!> \param nrow global row count
!> \param ncol global column count
!> \param fmat matrix read from restart
!> \param rst_unit restart unit
! **************************************************************************************************
   SUBROUTINE read_restart_fm_matrix(fmat_template, nrow, ncol, fmat, rst_unit)
      TYPE(cp_fm_type), INTENT(IN)                       :: fmat_template
      INTEGER, INTENT(IN)                                :: nrow, ncol
      TYPE(cp_fm_type), POINTER                          :: fmat
      INTEGER, INTENT(IN)                                :: rst_unit

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      NULLIFY (fm_struct)
      ALLOCATE (fmat)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=fmat_template%matrix_struct, &
                               nrow_global=nrow, ncol_global=ncol)
      CALL cp_fm_create(fmat, fm_struct)
      CALL cp_fm_read_unformatted(fmat, rst_unit)
      CALL cp_fm_struct_release(fm_struct)

   END SUBROUTINE read_restart_fm_matrix

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param cell_to_index ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_index(cell, cell_to_index) RESULT(index)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER                                            :: index

      IF (CELL(1) >= LBOUND(cell_to_index, 1) .AND. CELL(1) <= UBOUND(cell_to_index, 1) .AND. &
          CELL(2) >= LBOUND(cell_to_index, 2) .AND. CELL(2) <= UBOUND(cell_to_index, 2) .AND. &
          CELL(3) >= LBOUND(cell_to_index, 3) .AND. CELL(3) <= UBOUND(cell_to_index, 3)) THEN

         index = cell_to_index(cell(1), cell(2), cell(3))

      ELSE

         index = 0

      END IF

   END FUNCTION get_index

END MODULE kpoint_io
