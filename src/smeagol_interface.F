!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K+SMEAGOL interface.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************
MODULE smeagol_interface
   USE bibliography, ONLY: Ahart2024, &
                           Bailey2006, &
                           cite_reference
   USE cell_types, ONLY: cell_type, &
                         real_to_scaled, &
                         scaled_to_real
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_operations, ONLY: dbcsr_allocate_matrix_set
   USE cp_files, ONLY: close_file, &
                       open_file
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE cp_output_handling, ONLY: cp_get_iter_level_by_name, &
                                 cp_get_iter_nr
   USE cp_dbcsr_api, ONLY: dbcsr_copy, &
                           dbcsr_p_type, &
                           dbcsr_set
   USE input_constants, ONLY: smeagol_bulklead_left, &
                              smeagol_bulklead_leftright, &
                              smeagol_bulklead_right, &
                              smeagol_runtype_bulktransport, &
                              smeagol_runtype_emtransport
   USE kahan_sum, ONLY: accurate_dot_product, &
                        accurate_sum
   USE kinds, ONLY: dp, &
                    int_8
   USE kpoint_types, ONLY: get_kpoint_info, &
                           kpoint_type
   USE message_passing, ONLY: mp_para_env_type
#if defined(__SMEAGOL)
   USE mmpi_negf, ONLY: create_communicators_negf, &
                        destroy_communicators_negf
   USE mnegf_interface, ONLY: negf_interface
   USE negfmod, ONLY: smeagolglobal_em_nas => em_nas, &
                      smeagolglobal_em_nau => em_nau, &
                      smeagolglobal_em_nso => em_nso, &
                      smeagolglobal_em_nuo => em_nuo, &
                      smeagolglobal_negfon => negfon
#endif
   USE physcon, ONLY: bohr, &
                      evolt
   USE pw_grid_types, ONLY: pw_grid_type
   USE pw_types, ONLY: pw_r3d_rs_type
   USE qs_matrix_w, ONLY: compute_matrix_w
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_ks_types, ONLY: qs_ks_env_type, &
                          set_ks_env
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
   USE qs_rho_types, ONLY: qs_rho_get, &
                           qs_rho_type
   USE qs_subsys_types, ONLY: qs_subsys_get, &
                              qs_subsys_type
   USE scf_control_types, ONLY: scf_control_type
   USE smeagol_control_types, ONLY: smeagol_control_type
   USE smeagol_emtoptions, ONLY: ReadOptionsNEGF_DFT, &
                                 emtrans_deallocate_global_arrays, &
                                 emtrans_options, &
                                 reademtr
   USE smeagol_matrix_utils, ONLY: convert_dbcsr_to_distributed_siesta, &
                                   convert_distributed_siesta_to_dbcsr, &
                                   siesta_distrib_csc_struct_type, &
                                   siesta_struct_create, &
                                   siesta_struct_release
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smeagol_interface'

   PUBLIC :: run_smeagol_bulktrans, run_smeagol_emtrans
   PUBLIC :: smeagol_shift_v_hartree

CONTAINS

! **************************************************************************************************
!> \brief Save overlap, Kohn-Sham, electron density, and energy-density matrices of semi-infinite
!>        electrodes in SIESTA format.
!> \param qs_env  QuickStep environment
! **************************************************************************************************
   SUBROUTINE run_smeagol_bulktrans(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'run_smeagol_bulktrans'

      CHARACTER(len=32)                                  :: hst_fmt
      INTEGER                                            :: funit, handle, img, ispin, lead_label, &
                                                            log_unit, nimages, nspin
      INTEGER(kind=int_8)                                :: ielem
      INTEGER, DIMENSION(2)                              :: max_ij_cell_image
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints, has_unit_metric, not_regtest
      REAL(kind=dp)                                      :: H_to_Ry
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: matrix_siesta_1d, nelectrons
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: matrix_siesta_2d
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: matrix_kp_generic
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp, matrix_w_kp, &
                                                            rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(pw_r3d_rs_type), POINTER                      :: v_hartree_rspace
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(siesta_distrib_csc_struct_type)               :: siesta_struct
      TYPE(smeagol_control_type), POINTER                :: smeagol_control

      CALL get_qs_env(qs_env, dft_control=dft_control)
      smeagol_control => dft_control%smeagol_control
      IF (.NOT. (smeagol_control%smeagol_enabled .AND. smeagol_control%run_type == smeagol_runtype_bulktransport)) RETURN

      CALL timeset(routineN, handle)

      log_unit = cp_logger_get_default_io_unit()
      H_to_Ry = smeagol_control%to_smeagol_energy_units
      not_regtest = .NOT. smeagol_control%do_regtest

      lead_label = smeagol_control%lead_label
      nspin = dft_control%nspins

      NULLIFY (v_hartree_rspace)
      CALL get_qs_env(qs_env, energy=energy, para_env=para_env, subsys=subsys, &
                      scf_control=scf_control, &
                      do_kpoints=do_kpoints, kpoints=kpoints, &
                      matrix_ks_kp=matrix_ks_kp, matrix_s_kp=matrix_s_kp, &
                      rho=rho, sab_orb=sab_nl, v_hartree_rspace=v_hartree_rspace)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      IF (not_regtest) THEN
         ! save average electrostatic potential of the electrode along transport direction
         CALL write_average_hartree_potential(v_hartree_rspace, smeagol_control%project_name)
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,T61,ES20.10E2)') "SMEAGOL| E_FERMI [a.u.] = ", energy%efermi
      END IF

      IF (do_kpoints) THEN
         nimages = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimages = 1
         ALLOCATE (cell_to_index(0:0, 0:0, 0:0))
         cell_to_index(0, 0, 0) = 1
         ! We do need at least two cell images along transport direction.
         CPABORT("Please enable k-points")
      END IF

      ! largest index of cell images along i and j cell vectors
      ! e.g., (2,0) in case of 5 cell images (0,0), (1,0), (-1,0), (2,0), and (-2,0)
      ! -1 means to use all available cell images along a particular cell vector.
      max_ij_cell_image(:) = -1
      DO img = 1, 2
         IF (smeagol_control%n_cell_images(img) > 0) THEN
            max_ij_cell_image(img) = smeagol_control%n_cell_images(img)/2
         END IF
      END DO

      ALLOCATE (matrix_kp_generic(nimages))

      ! compute energy-density (W) matrix. We may need it later to calculate NEGF forces
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)
      IF (.NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, matrix_w_kp=matrix_w_kp)
         IF (.NOT. ASSOCIATED(matrix_w_kp)) THEN
            NULLIFY (matrix_w_kp)
            CALL get_qs_env(qs_env, ks_env=ks_env)
            CALL dbcsr_allocate_matrix_set(matrix_w_kp, nspin, nimages)
            DO ispin = 1, nspin
               DO img = 1, nimages
                  ALLOCATE (matrix_w_kp(ispin, img)%matrix)
                  CALL dbcsr_copy(matrix_w_kp(ispin, img)%matrix, matrix_s_kp(1, 1)%matrix, name="W MATRIX")
                  CALL dbcsr_set(matrix_w_kp(ispin, img)%matrix, 0.0_dp)
               END DO
            END DO
            CALL set_ks_env(ks_env, matrix_w_kp=matrix_w_kp)
         END IF
      END IF
      CALL compute_matrix_w(qs_env, calc_forces=.TRUE.)

      ! obtain the sparsity pattern of the overlap matrix
      DO img = 1, nimages
         matrix_kp_generic(img)%matrix => matrix_s_kp(1, img)%matrix
      END DO

      CALL siesta_struct_create(siesta_struct, matrix_kp_generic, subsys, cell_to_index, &
                                sab_nl, para_env, max_ij_cell_image, do_merge=.FALSE., gather_root=para_env%source)

      ! write 'bulklft.DAT' and 'bulkrgt.DAT' files
      funit = -1
      IF (not_regtest .AND. para_env%mepos == para_env%source) THEN
         ! I/O process
         IF (lead_label == smeagol_bulklead_left .OR. lead_label == smeagol_bulklead_leftright) THEN
            CALL open_file("bulklft.DAT", file_status="REPLACE", file_form="FORMATTED", file_action="WRITE", unit_number=funit)

            CALL write_bulk_dat_file(funit, siesta_struct, smeagol_control%project_name, nspin, &
                                     energy%efermi, scf_control%smear%ELECTRONIC_TEMPERATURE, &
                                     H_to_Ry, do_kpoints, max_ij_cell_image)

            CALL close_file(funit)
         END IF

         IF (lead_label == smeagol_bulklead_right .OR. lead_label == smeagol_bulklead_leftright) THEN
            CALL open_file("bulkrgt.DAT", file_status="REPLACE", file_form="FORMATTED", file_action="WRITE", unit_number=funit)

            CALL write_bulk_dat_file(funit, siesta_struct, smeagol_control%project_name, nspin, &
                                     energy%efermi, scf_control%smear%ELECTRONIC_TEMPERATURE, &
                                     H_to_Ry, do_kpoints, max_ij_cell_image)

            CALL close_file(funit)
         END IF
      END IF

      ! write project_name.HST file
      funit = -1
      IF (para_env%mepos == para_env%source) THEN
         ALLOCATE (matrix_siesta_1d(siesta_struct%n_nonzero_elements))
         ALLOCATE (matrix_siesta_2d(siesta_struct%n_nonzero_elements, nspin))

         IF (not_regtest) THEN
            CALL open_file(TRIM(smeagol_control%project_name)//".HST", &
                           file_status="REPLACE", file_form="FORMATTED", file_action="WRITE", unit_number=funit)
         END IF
      ELSE
         ALLOCATE (matrix_siesta_1d(1))
         ALLOCATE (matrix_siesta_2d(1, nspin))
      END IF

      !DO img = 1, nimages
      !   matrix_kp_generic(img)%matrix => matrix_s_kp(1, img)%matrix
      !END DO
      CALL convert_dbcsr_to_distributed_siesta(matrix_siesta_1d, matrix_kp_generic, siesta_struct, para_env)

      DO ispin = 1, nspin
         DO img = 1, nimages
            matrix_kp_generic(img)%matrix => matrix_ks_kp(ispin, img)%matrix
         END DO
         CALL convert_dbcsr_to_distributed_siesta(matrix_siesta_2d(:, ispin), matrix_kp_generic, siesta_struct, para_env)
      END DO
      ! As SIESTA's default energy unit is Rydberg, scale the KS-matrix
      matrix_siesta_2d(:, :) = H_to_Ry*matrix_siesta_2d(:, :)

      IF (funit > 0) THEN ! not_regtest .AND. para_env%mepos == para_env%source
         WRITE (hst_fmt, '(A,I0,A)') "(", nspin, "ES26.17E3)"
         DO ielem = 1, siesta_struct%n_nonzero_elements
            WRITE (funit, '(ES26.17E3)') matrix_siesta_1d(ielem)
            WRITE (funit, hst_fmt) matrix_siesta_2d(ielem, :)
         END DO

         CALL close_file(funit)
      END IF

      ! write density matrix
      DO ispin = 1, nspin
         DO img = 1, nimages
            matrix_kp_generic(img)%matrix => rho_ao_kp(ispin, img)%matrix
         END DO
         CALL convert_dbcsr_to_distributed_siesta(matrix_siesta_2d(:, ispin), matrix_kp_generic, siesta_struct, para_env)
      END DO

      IF (para_env%mepos == para_env%source) THEN
         ALLOCATE (nelectrons(nspin))
         DO ispin = 1, nspin
            nelectrons(ispin) = accurate_dot_product(matrix_siesta_2d(:, ispin), matrix_siesta_1d)
         END DO

         CPASSERT(log_unit > 0)
         IF (nspin > 1) THEN
            WRITE (log_unit, '(T2,A,T61,F20.10)') "SMEAGOL| Number of alpha-spin electrons: ", nelectrons(1)
            WRITE (log_unit, '(T2,A,T61,F20.10)') "SMEAGOL| Number of  beta-spin electrons: ", nelectrons(2)
         ELSE
            WRITE (log_unit, '(T2,A,T61,F20.10)') "SMEAGOL| Number of electrons: ", nelectrons(1)
         END IF
         DEALLOCATE (nelectrons)

         IF (not_regtest) THEN
            CALL open_file(TRIM(smeagol_control%project_name)//".DM", &
                           file_status="REPLACE", file_form="UNFORMATTED", file_action="WRITE", unit_number=funit)

            CALL write_bulk_dm_file(funit, siesta_struct, nspin, matrix_siesta_2d)

            CALL close_file(funit)
         END IF
      END IF

      CALL get_qs_env(qs_env, matrix_w_kp=matrix_w_kp)
      IF (ASSOCIATED(matrix_w_kp)) THEN
         ! write energy density matrix
         DO ispin = 1, nspin
            DO img = 1, nimages
               matrix_kp_generic(img)%matrix => matrix_w_kp(ispin, img)%matrix
            END DO
            CALL convert_dbcsr_to_distributed_siesta(matrix_siesta_2d(:, ispin), matrix_kp_generic, &
                                                     siesta_struct, para_env)
         END DO

         IF (not_regtest .AND. para_env%mepos == para_env%source) THEN
            CALL open_file(TRIM(smeagol_control%project_name)//".EDM", &
                           file_status="REPLACE", file_form="UNFORMATTED", file_action="WRITE", unit_number=funit)

            CALL write_bulk_dm_file(funit, siesta_struct, nspin, matrix_siesta_2d)

            CALL close_file(funit)
         END IF
      END IF

      DEALLOCATE (matrix_siesta_2d)
      DEALLOCATE (matrix_siesta_1d)

      DEALLOCATE (matrix_kp_generic)
      IF (.NOT. do_kpoints) DEALLOCATE (cell_to_index)

      CALL siesta_struct_release(siesta_struct)
      CALL timestop(handle)
   END SUBROUTINE run_smeagol_bulktrans

! **************************************************************************************************
!> \brief Write a sparse matrix structure file in SIESTA format. Should be called on I/O MPI process only.
!> \param funit             file to write
!> \param siesta_struct     sparse matrix structure in SIESTA format
!> \param system_label      SMEAGOL project name (first components of file names, i.e. system_label.HST)
!> \param nspin             number of spin components
!> \param EFermi            Fermi level
!> \param temperature       electronic temperature
!> \param H_to_Ry           Hartree to Rydberg scale factor
!> \param do_kpoints        whether to perform k-point calculation. Should always be enabled as
!>                          SMEAGOL expects at least 3 cell replicas along the transport direction
!> \param max_ij_cell_image maximum cell-replica indices along i- and j- lattice vectors
!>                          (perpendicular the transport direction)
! **************************************************************************************************
   SUBROUTINE write_bulk_dat_file(funit, siesta_struct, system_label, nspin, EFermi, temperature, &
                                  H_to_Ry, do_kpoints, max_ij_cell_image)
      INTEGER, INTENT(in)                                :: funit
      TYPE(siesta_distrib_csc_struct_type), INTENT(in)   :: siesta_struct
      CHARACTER(len=*), INTENT(in)                       :: system_label
      INTEGER, INTENT(in)                                :: nspin
      REAL(kind=dp), INTENT(in)                          :: EFermi, temperature, H_to_Ry
      LOGICAL, INTENT(in)                                :: do_kpoints
      INTEGER, DIMENSION(2), INTENT(in)                  :: max_ij_cell_image

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_bulk_dat_file'

      INTEGER                                            :: handle, icol, irow, nao_supercell, &
                                                            nao_unitcell
      INTEGER, DIMENSION(2)                              :: ncells_siesta

      CALL timeset(routineN, handle)

      ! ++ header
      nao_unitcell = siesta_struct%nrows
      nao_supercell = siesta_struct%ncols
      ncells_siesta(1:2) = 2*max_ij_cell_image(1:2) + 1

      ! SMEAGOL expects Temperature and Fermi energy in Rydberg energy units, not in Hartree energy units.
      ! It is why these values are doubled.

      WRITE (funit, '(1X,A20,3I12,2ES26.17E3,3I12)') &
         system_label, nao_unitcell, nspin, siesta_struct%n_nonzero_elements, &
         H_to_Ry*EFermi, H_to_Ry*temperature, ncells_siesta(1:2), nao_supercell

      ! ++ number of non-zero matrix elements on each row and
      !    the index of the first non-zero matrix element on this row -1 in 1D data array.
      DO irow = 1, nao_unitcell
         WRITE (funit, '(2I12)') siesta_struct%n_nonzero_cols(irow), siesta_struct%row_offset(irow)
      END DO

      DO icol = 1, nao_supercell
         WRITE (funit, '(I12)') siesta_struct%indxuo(icol)
      END DO

      ! ++ column indices of non-zero matrix elements
      DO irow = 1, nao_unitcell
         DO icol = 1, siesta_struct%n_nonzero_cols(irow)
            WRITE (funit, '(I12)') siesta_struct%col_index(siesta_struct%row_offset(irow) + icol)

            IF (do_kpoints) THEN
               WRITE (funit, '(F21.16,5X,F21.16,5X,F21.16)') siesta_struct%xij(:, siesta_struct%row_offset(irow) + icol)
            END IF
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE write_bulk_dat_file

! **************************************************************************************************
!> \brief Write an (energy-) density matrix. Should be called on I/O MPI process only.
!> \param funit          file to write
!> \param siesta_struct  sparse matrix structure in SIESTA format
!> \param nspin          number of spin components
!> \param matrix_siesta  non-zero matrix elements (1:siesta_struct%n_nonzero_elements, 1:nspin)
! **************************************************************************************************
   SUBROUTINE write_bulk_dm_file(funit, siesta_struct, nspin, matrix_siesta)
      INTEGER, INTENT(in)                                :: funit
      TYPE(siesta_distrib_csc_struct_type), INTENT(in)   :: siesta_struct
      INTEGER, INTENT(in)                                :: nspin
      REAL(kind=dp), DIMENSION(:, :), INTENT(in)         :: matrix_siesta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_bulk_dm_file'

      INTEGER                                            :: handle, irow, ispin

      CALL timeset(routineN, handle)

      ! ++ number of compressed rows, number of spin components
      WRITE (funit) siesta_struct%nrows, nspin

      ! ++ number of non-zero matrix elements on each compressed row.
      !    The sparsity pattern for alpha- and beta-spin density matrices are identical
      WRITE (funit) siesta_struct%n_nonzero_cols

      ! ++ column indices of non-zero matrix elements
      DO irow = 1, siesta_struct%nrows
         WRITE (funit) siesta_struct%col_index(siesta_struct%row_offset(irow) + 1: &
                                               siesta_struct%row_offset(irow) + siesta_struct%n_nonzero_cols(irow))
      END DO

      ! ++ non-zero matrix blocks
      DO ispin = 1, nspin
         DO irow = 1, siesta_struct%nrows
            WRITE (funit) matrix_siesta(siesta_struct%row_offset(irow) + 1: &
                                        siesta_struct%row_offset(irow) + siesta_struct%n_nonzero_cols(irow), ispin)
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE write_bulk_dm_file

! **************************************************************************************************
!> \brief Write the average value of Hartree potential along transport direction.
!>        SMEAGOL assumes that the transport direction coincides with z-axis.
!> \param v_hartree_rspace   Hartree potential on a real-space 3-D grid
!> \param project_name       SMEAGOL project name
!> \note This routine assumes that the lattice vector 'k' coincides with z-axis
! **************************************************************************************************
   SUBROUTINE write_average_hartree_potential(v_hartree_rspace, project_name)
      TYPE(pw_r3d_rs_type), POINTER                      :: v_hartree_rspace
      CHARACTER(len=*), INTENT(in)                       :: project_name

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_average_hartree_potential'

      INTEGER                                            :: funit, handle, iz, lz, uz
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: v_hartree_z_average
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: cr3d
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CALL timeset(routineN, handle)

      pw_grid => v_hartree_rspace%pw_grid
      cr3d => v_hartree_rspace%array

      ALLOCATE (v_hartree_z_average(pw_grid%bounds(1, 3):pw_grid%bounds(2, 3)))
      v_hartree_z_average(:) = 0.0_dp

      lz = pw_grid%bounds_local(1, 3)
      uz = pw_grid%bounds_local(2, 3)

      ! save average electrostatic potential
      DO iz = lz, uz
         v_hartree_z_average(iz) = accurate_sum(cr3d(:, :, iz))
      END DO
      CALL pw_grid%para%group%sum(v_hartree_z_average)
      v_hartree_z_average(:) = v_hartree_z_average(:)/ &
                               (REAL(pw_grid%npts(1), kind=dp)*REAL(pw_grid%npts(2), kind=dp))

      funit = -1
      IF (pw_grid%para%group%mepos == pw_grid%para%group%source) THEN
         CALL open_file(TRIM(ADJUSTL(project_name))//"-VH_AV.dat", &
                        file_status="REPLACE", file_form="FORMATTED", file_action="WRITE", unit_number=funit)
         WRITE (funit, '(A,T10,A,T25,A)') "#", "z (A)", "V_H average (eV)"
         DO iz = lz, uz
            WRITE (funit, '(F20.10,ES20.10E3)') pw_grid%dh(3, 3)*REAL(iz - lz, kind=dp)/bohr, &
               v_hartree_z_average(iz)/pw_grid%dvol*evolt
         END DO
         CALL close_file(funit)
      END IF

      DEALLOCATE (v_hartree_z_average)
      CALL timestop(handle)
   END SUBROUTINE write_average_hartree_potential

! **************************************************************************************************
!> \brief Align Hatree potential of semi-infinite leads to match bulk-transport calculation
!>        and apply external electrostatic potential (bias).
!> \param v_hartree_rspace     Hartree potential on a real-space 3-D grid [inout]
!> \param cell                 unit cell
!> \param HartreeLeadsLeft     z-coordinate of the left lead
!> \param HartreeLeadsRight    z-coordinate of the right lead
!> \param HartreeLeadsBottom   average Hartree potential (from bulk-transport calculation)
!>                             at HartreeLeadsLeft and HartreeLeadsRight
!> \param Vbias                the value of external potential to apply
!> \param zleft                starting point of external potential drop (initial value 0.5*Vbias)
!> \param zright               final point of external potential drop (final value -0.5*Vbias)
!> \param isexplicit_zright    whether zright has beed provided explicitly via the input file.
!>                             If not, use the cell boundary.
!> \param isexplicit_bottom    whether the reference Hatree potential for bulk regions has been
!>                             provided explicitly via the input file. If not, do not align the
!>                             potential at all (instead of aligning it to 0 which is incorrect).
!> \note This routine assumes that the lattice vector 'k' coincides with z-axis
! **************************************************************************************************
   SUBROUTINE smeagol_shift_v_hartree(v_hartree_rspace, cell, &
                                      HartreeLeadsLeft, HartreeLeadsRight, HartreeLeadsBottom, &
                                      Vbias, zleft, zright, isexplicit_zright, isexplicit_bottom)
      TYPE(pw_r3d_rs_type), POINTER                      :: v_hartree_rspace
      TYPE(cell_type), POINTER                           :: cell
      REAL(kind=dp), INTENT(in)                          :: HartreeLeadsLeft, HartreeLeadsRight, &
                                                            HartreeLeadsBottom, Vbias, zleft, &
                                                            zright
      LOGICAL, INTENT(in)                                :: isexplicit_zright, isexplicit_bottom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'smeagol_shift_v_hartree'

      INTEGER                                            :: handle, iz, l_right, lz, u_left, uz
      REAL(kind=dp)                                      :: v_average_left, v_average_right, &
                                                            v_bias_iz, zright_explicit
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: v_hartree_z_average
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: cr3d
      REAL(kind=dp), DIMENSION(3)                        :: r, r_pbc
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CALL timeset(routineN, handle)
      pw_grid => v_hartree_rspace%pw_grid
      cr3d => v_hartree_rspace%array

      zright_explicit = zright
      IF (.NOT. isexplicit_zright) THEN
         r_pbc(:) = (/0.0_dp, 0.0_dp, 1.0_dp/)
         CALL scaled_to_real(r, r_pbc, cell)
         zright_explicit = r(3)
      END IF

      ALLOCATE (v_hartree_z_average(pw_grid%bounds(1, 3):pw_grid%bounds(2, 3)))

      lz = pw_grid%bounds_local(1, 3)
      uz = pw_grid%bounds_local(2, 3)

      v_hartree_z_average(:) = 0.0_dp
      DO iz = lz, uz
         v_hartree_z_average(iz) = accurate_sum(cr3d(:, :, iz))
      END DO

      CALL pw_grid%para%group%sum(v_hartree_z_average)
      v_hartree_z_average(:) = v_hartree_z_average(:)/ &
                               (REAL(pw_grid%npts(1), kind=dp)*REAL(pw_grid%npts(2), kind=dp))

      ! z-indices of the V_hartree related to the left lead: pw_grid%bounds(1,3) .. u_left
      r(1:3) = (/0.0_dp, 0.0_dp, HartreeLeadsLeft/)
      CALL real_to_scaled(r_pbc, r, cell)
      u_left = NINT(r_pbc(3)*REAL(pw_grid%npts(3), kind=dp)) + pw_grid%bounds(1, 3)
      IF (u_left > pw_grid%bounds(2, 3)) u_left = pw_grid%bounds(2, 3)

      ! z-indices of the V_hartree related to the right lead: l_right .. pw_grid%bounds(2, 3)
      r(1:3) = (/0.0_dp, 0.0_dp, HartreeLeadsRight/)
      CALL real_to_scaled(r_pbc, r, cell)
      l_right = NINT(r_pbc(3)*REAL(pw_grid%npts(3), kind=dp)) + pw_grid%bounds(1, 3)
      IF (l_right > pw_grid%bounds(2, 3)) l_right = pw_grid%bounds(2, 3)

      CPASSERT(u_left <= l_right)

      v_average_left = v_hartree_z_average(u_left)
      v_average_right = v_hartree_z_average(l_right)

      ! align electrostatic potential of leads' regions with ones from bulk transport calculation
      IF (isexplicit_bottom) THEN
         v_hartree_z_average(:) = HartreeLeadsBottom*pw_grid%dvol - 0.5_dp*(v_average_left + v_average_right)
      ELSE
         ! do not align electrostatic potential
         v_hartree_z_average(:) = 0.0_dp
      END IF

      ! external Vbias
      ! TO DO: convert zright and zleft to scaled coordinates instead
      DO iz = lz, uz
         r_pbc(1) = 0.0_dp
         r_pbc(2) = 0.0_dp
         r_pbc(3) = REAL(iz - pw_grid%bounds(1, 3), kind=dp)/REAL(pw_grid%npts(3), kind=dp)
         CALL scaled_to_real(r, r_pbc, cell)
         IF (r(3) < zleft) THEN
            v_bias_iz = 0.5_dp*Vbias
         ELSE IF (r(3) > zright_explicit) THEN
            v_bias_iz = -0.5_dp*Vbias
         ELSE
            v_bias_iz = Vbias*(0.5_dp - (r(3) - zleft)/(zright_explicit - zleft))
         END IF
         v_hartree_z_average(iz) = v_hartree_z_average(iz) + v_bias_iz*pw_grid%dvol
      END DO

      DO iz = lz, uz
         cr3d(:, :, iz) = cr3d(:, :, iz) + v_hartree_z_average(iz)
      END DO

      DEALLOCATE (v_hartree_z_average)
      CALL timestop(handle)
   END SUBROUTINE smeagol_shift_v_hartree

! **************************************************************************************************
!> \brief Run NEGF/SMEAGOL transport calculation.
!> \param qs_env     QuickStep environment
!> \param last       converged SCF iterations; compute NEGF properties [in]
!> \param iter       index of the current iteration [in]
!> \param rho_ao_kp  refined electron density; to be mixed with electron density from the previous
!>                   SCF iteration [out]
! **************************************************************************************************
   SUBROUTINE run_smeagol_emtrans(qs_env, last, iter, rho_ao_kp)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: last
      INTEGER, INTENT(in)                                :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: rho_ao_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'run_smeagol_emtrans'

      INTEGER                                            :: handle, img, ispin, md_step, natoms, &
                                                            nimages, nspin
      INTEGER, DIMENSION(2)                              :: max_ij_cell_image, n_cell_images
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints, local_ldos, local_TrCoeff, &
                                                            negfon_saved, structure_changed
      REAL(kind=dp)                                      :: H_to_Ry
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: matrix_s_csc_merged
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: Dnew_csc_merged, Enew_csc_merged, &
                                                            matrix_ks_csc_merged
      TYPE(cell_type), POINTER                           :: ucell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: matrix_kp_generic
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp, matrix_w_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(siesta_distrib_csc_struct_type)               :: siesta_struct_merged
      TYPE(smeagol_control_type), POINTER                :: smeagol_control

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, dft_control=dft_control)

      nspin = dft_control%nspins
      smeagol_control => dft_control%smeagol_control
      H_to_Ry = smeagol_control%to_smeagol_energy_units

      IF (.NOT. (smeagol_control%smeagol_enabled .AND. smeagol_control%run_type == smeagol_runtype_emtransport)) RETURN
      CALL timeset(routineN, handle)

      NULLIFY (kpoints, matrix_s_kp, matrix_ks_kp, matrix_w_kp, para_env, sab_nl, scf_control, subsys)
#if defined(__SMEAGOL)
      CALL cite_reference(Ahart2024)
      CALL cite_reference(Bailey2006)

      CPASSERT(ASSOCIATED(smeagol_control%aux))
      CALL get_qs_env(qs_env, para_env=para_env, scf_control=scf_control, &
                      do_kpoints=do_kpoints, kpoints=kpoints, &
                      matrix_ks_kp=matrix_ks_kp, matrix_s_kp=matrix_s_kp, matrix_w_kp=matrix_w_kp, &
                      sab_orb=sab_nl, subsys=subsys)
      CALL qs_subsys_get(subsys, cell=ucell, natom=natoms)

      IF (do_kpoints) THEN
         nimages = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimages = 1
         ALLOCATE (cell_to_index(0:0, 0:0, 0:0))
         cell_to_index(0, 0, 0) = 1
      END IF

      max_ij_cell_image(:) = -1
      DO img = 1, 2
         IF (smeagol_control%n_cell_images(img) > 0) THEN
            max_ij_cell_image(img) = smeagol_control%n_cell_images(img)/2
         END IF
      END DO

      ! obtain the sparsity pattern of the Kohn-Sham matrix
      ALLOCATE (matrix_kp_generic(nimages))
      DO img = 1, nimages
         matrix_kp_generic(img)%matrix => matrix_s_kp(1, img)%matrix
      END DO

      CALL siesta_struct_create(siesta_struct_merged, matrix_kp_generic, subsys, &
                                cell_to_index, sab_nl, para_env, max_ij_cell_image, do_merge=.TRUE., gather_root=-1)

      ! Number of unit cells along (x and y ?) directions
      n_cell_images(1:2) = 2*max_ij_cell_image(1:2) + 1

      ALLOCATE (matrix_s_csc_merged(siesta_struct_merged%n_nonzero_elements))
      ALLOCATE (matrix_ks_csc_merged(siesta_struct_merged%n_nonzero_elements, nspin))

      CALL convert_dbcsr_to_distributed_siesta(matrix_s_csc_merged, matrix_kp_generic, siesta_struct_merged, para_env)

      DO ispin = 1, nspin
         DO img = 1, nimages
            matrix_kp_generic(img)%matrix => matrix_ks_kp(ispin, img)%matrix
         END DO

         CALL convert_dbcsr_to_distributed_siesta(matrix_ks_csc_merged(:, ispin), matrix_kp_generic, &
                                                  siesta_struct_merged, para_env)
      END DO

      IF (smeagol_control%aux%md_iter_level > 0) THEN
         CALL cp_get_iter_nr(logger%iter_info, smeagol_control%aux%md_iter_level, iter_nr=md_step)
      ELSE IF (smeagol_control%aux%md_iter_level == 0) THEN
         ! single-point energy calculation : there is only one MD step in this case
         md_step = smeagol_control%aux%md_first_step
      ELSE
         ! first invocation of the subroutine : initialise md_iter_level and md_first_step variables
         smeagol_control%aux%md_iter_level = cp_get_iter_level_by_name(logger%iter_info, "MD")

         IF (smeagol_control%aux%md_iter_level <= 0) &
            smeagol_control%aux%md_iter_level = cp_get_iter_level_by_name(logger%iter_info, "GEO_OPT")

         IF (smeagol_control%aux%md_iter_level <= 0) &
            smeagol_control%aux%md_iter_level = 0

         !  index of the first GEO_OPT / MD step
         IF (smeagol_control%aux%md_iter_level > 0) THEN
            CALL cp_get_iter_nr(logger%iter_info, smeagol_control%aux%md_iter_level, iter_nr=smeagol_control%aux%md_first_step)
         ELSE
            ! it has already been initialised in read_smeagol_control()
            smeagol_control%aux%md_first_step = 0
         END IF

         md_step = smeagol_control%aux%md_first_step
      END IF

      CALL reademtr(smeagol_control, natoms, gamma_negf=.FALSE.)
      CALL ReadOptionsNEGF_DFT(smeagol_control, ucell, torqueflag=.FALSE., torquelin=.FALSE.)

      CALL emtrans_options(smeagol_control, &
                           matrix_s=matrix_s_kp(1, 1)%matrix, para_env=para_env, iter=iter, &
                           istep=md_step, inicoor=smeagol_control%aux%md_first_step, iv=0, &
                           delta=smeagol_control%aux%delta, nk=kpoints%nkp)

      CALL create_communicators_negf(parent_comm=para_env%get_handle(), &
                                     nprocs_inverse=smeagol_control%aux%nprocs_inverse, &
                                     NParallelK=smeagol_control%aux%NParallelK)

      ALLOCATE (Dnew_csc_merged(siesta_struct_merged%n_nonzero_elements, nspin))
      ALLOCATE (Enew_csc_merged(siesta_struct_merged%n_nonzero_elements, nspin))

      ! As SMEAGOL's default energy unit is Rydberg, scale the KS-matrix
      matrix_ks_csc_merged(:, :) = H_to_Ry*matrix_ks_csc_merged(:, :)

      ! SMEAGOL computes current if EM.OrderN = .false., otherwise the printed current is equal to 0.
      ! The following computes current regardless the value of EM.OrderN parameter
      IF (last) THEN
         negfon_saved = smeagolglobal_negfon
         smeagolglobal_negfon = .FALSE.
      END IF

      IF (last) THEN
         local_TrCoeff = smeagol_control%aux%TrCoeff
         local_ldos = .TRUE. ! TO DO: find out initial value of ldos
      ELSE
         local_TrCoeff = .FALSE.
         local_ldos = .FALSE.
      END IF

      ! number of atoms in the unit cell
      smeagolglobal_em_nau = natoms ! number of atoms in the unit cell

      ! number of atomic orbitals in the unit cell.
      ! This global variable defines the size of temporary arrays for (P)DOS calculation.
      ! This should be the total number of the atomic orbitals in the unit cell,
      ! instead of the number of atomic orbitals local to the current MPI process.
      smeagolglobal_em_nuo = siesta_struct_merged%ncols/(n_cell_images(1)*n_cell_images(2))

      ! number of atoms in the super cell
      smeagolglobal_em_nas = SIZE(siesta_struct_merged%xa)

      ! number of atomic orbitals in the super cell
      smeagolglobal_em_nso = siesta_struct_merged%ncols

      ! The following global variables are only used in writephibin() which is not called by this interface
      !   smeagolglobal_em_isa => isa
      !   smeagolglobal_em_iaorb => iaorb
      !   smeagolglobal_em_iphorb => iphorb

      structure_changed = (iter == 1)

      CALL Negf_Interface( &
         ! distributed Kohn-Sham, overlap, density, and energy-density matrices in SIESTA format
         H=matrix_ks_csc_merged, &
         S=matrix_s_csc_merged, &
         DM=Dnew_csc_merged, &
         Omega=Enew_csc_merged, &
         ! interatomic distances for each non-zero matrix element
         xij=siesta_struct_merged%xij, &
         ! number of atomic orbitals in a supercell
         no_s=siesta_struct_merged%ncols, &
         ! number of atomic orbitals in the unit cell
         no_u=siesta_struct_merged%ncols/(n_cell_images(1)*n_cell_images(2)), &
         ! number of AOs local to the given MPI process
         no_u_node=siesta_struct_merged%nrows, &
         ! atomic coordinates for each atom in the supercell
         xa=siesta_struct_merged%xa, &
         ! unused dummy variable na_u
         na_u=natoms, &
         ! number of atoms in the supercell
         na_s=SIZE(siesta_struct_merged%xa), &
         ! number of spin-components
         NspinRealInputMatrix=nspin, &
         ! number of non-zero matrix element
         maxnh=INT(siesta_struct_merged%n_nonzero_elements), &
         ! number of non-zero matrix elements on each locally stored row
         numh=siesta_struct_merged%n_nonzero_cols, &
         ! offset (index-1) of first non-zero matrix elements on each locally stored row
         listhptr=siesta_struct_merged%row_offset, &
         ! column number (AO index) for each non-zero matrix element
         listh=siesta_struct_merged%col_index, &
         ! number of k-points
         nkpts=kpoints%nkp, &
         ! k-point coordinates
         kpoint=kpoints%xkp, &
         ! k-point weight
         weight_k=kpoints%wkp, &
         ! index of equivalent atomic orbital within the primary unit cell for each AO in the supercell
         indxuo=siesta_struct_merged%indxuo, &
         ! list of atomic indices on which each AO (in the supercell) is centred
         iaorb=siesta_struct_merged%iaorb, &
         ! GEO_OPT / MD step index
         MD_Step=md_step, &
         ! GEO_OPT / MD at which SMEAGOL should allocate its internal data structures
         inicoor=smeagol_control%aux%md_first_step, &
         ! ivv (step over bias, first IV step should always be 0 (hardcoded in SMEAGOL)
         !    we iterate over GEO_OPT / MD steps instead of bias steps in order not to overwrite
         !    TRC (Transmission coefficients) files
         IV_step=md_step - smeagol_control%aux%md_first_step, &
         ! applied electrostatic bias
         Vb=smeagol_control%aux%VBias*H_to_Ry, &
         ! index of the currect SCF iteration
         SCF_step=iter, &
         ! compute density matrix (.FALSE.) or properties (.TRUE.)
         Last_SCF_step=last, &
         ! recompute self-energy matrices
         StructureChanged=structure_changed, &
         ! electronic temperature in Ry
         temp=smeagol_control%aux%temperature*H_to_Ry, &
         ! number of unit cell replicas along i-, and j- cell verctors
         nsc=n_cell_images, &
         ! name of SMEAGOL project (partial file name for files created by SMEAGOL)
         slabel=smeagol_control%project_name, &
         ! number of integration points along real axis, circular path and a line in imaginary space parallel to the real axis
         NEnergR=smeagol_control%aux%NEnergR, &
         NEnergIC=smeagol_control%aux%NEnergIC, &
         NEnergIL=smeagol_control%aux%NEnergIL, &
         ! number of poles
         NPoles=smeagol_control%aux%NPoles, &
         ! small imaginary shift that makes matrix inversion computationally stable
         Delta=smeagol_control%aux%deltamin*H_to_Ry, &
         ! integration lower bound
         EnergLB=smeagol_control%aux%EnergLB*H_to_Ry, &
         ! [unused dummy arguments] initial (VInitial) and final (VFinal) voltage
         VInitial=smeagol_control%aux%VBias, &
         VFinal=smeagol_control%aux%VBias, &
         !
         SpinCL=smeagol_control%aux%SpinCL, &
         ! number of slices for OrderN matrix inversion
         NSlices=smeagol_control%aux%NSlices, &
         ! whether to compute transmission coefficients (TrCoeff), IETS spectrum (CalcIETS), and local Dnsity-of-States (ldos)
         TrCoeff=local_TrCoeff, &
         CalcIETS=smeagol_control%aux%CalcIETS, &
         ldos=local_ldos, &
         ! Transmission coefficients are only computed for certain runtypes (encoded with magic numbers).
         ! In case of idyn=0, transmission coefficients are always computed.
         idyn=0, &
         ! do not compute transmission coefficient for initial GEO_OPT / MD iterations
         tmdskip=smeagol_control%aux%tmdskip, &
         ! computes transmission coefficients once for each 'tmdsampling' MD iterations
         tmdsampling=smeagol_control%aux%tmdsampling)

      ! *** Bound-state correction method ***

      ! smeagol_control%bs_add    : BS.Add (.FALSE.) ; add bound states
      ! smeagol_control%bs_method : BS.Method (0) ; 0 - use effective Hamiltonian; 1 - add small imaginary part to the selfenergies
      ! smeagol_control%bssc      : BS.SetOccupation (1)
      IF (smeagol_control%aux%bs_add .AND. smeagol_control%aux%bs_method == 1 .AND. smeagol_control%aux%bssc == 0 .AND. &
          kpoints%nkp > 1 .AND. (.NOT. last)) THEN

         CALL Negf_Interface(H=matrix_ks_csc_merged, &
                             S=matrix_s_csc_merged, &
                             DM=Dnew_csc_merged, &
                             Omega=Enew_csc_merged, &
                             xij=siesta_struct_merged%xij, &
                             no_s=siesta_struct_merged%ncols, &
                             no_u=siesta_struct_merged%ncols/(n_cell_images(1)*n_cell_images(2)), &
                             no_u_node=siesta_struct_merged%nrows, &
                             xa=siesta_struct_merged%xa, &
                             ! unused dummy variable
                             na_u=natoms, &
                             na_s=SIZE(siesta_struct_merged%xa), &
                             NspinRealInputMatrix=nspin, &
                             maxnh=INT(siesta_struct_merged%n_nonzero_elements), &
                             numh=siesta_struct_merged%n_nonzero_cols, &
                             listhptr=siesta_struct_merged%row_offset, &
                             listh=siesta_struct_merged%col_index, &
                             nkpts=kpoints%nkp, &
                             kpoint=kpoints%xkp, &
                             weight_k=kpoints%wkp, &
                             indxuo=siesta_struct_merged%indxuo, &
                             iaorb=siesta_struct_merged%iaorb, &
                             MD_Step=md_step, &
                             inicoor=smeagol_control%aux%md_first_step, &
                             IV_step=md_step - smeagol_control%aux%md_first_step, &
                             Vb=smeagol_control%aux%VBias*H_to_Ry, &
                             ! This line is the only difference from the first Negf_Interface() call
                             SCF_step=MERGE(2, 1, structure_changed), &
                             Last_SCF_step=last, &
                             StructureChanged=structure_changed, &
                             temp=smeagol_control%aux%temperature*H_to_Ry, &
                             nsc=n_cell_images, &
                             slabel=smeagol_control%project_name, &
                             NEnergR=smeagol_control%aux%NEnergR, &
                             NEnergIC=smeagol_control%aux%NEnergIC, &
                             NEnergIL=smeagol_control%aux%NEnergIL, &
                             NPoles=smeagol_control%aux%NPoles, &
                             Delta=smeagol_control%aux%deltamin*H_to_Ry, &
                             EnergLB=smeagol_control%aux%EnergLB*H_to_Ry, &
                             ! unused dummy variable
                             VInitial=smeagol_control%aux%VBias, &
                             ! unused dummy variable
                             VFinal=smeagol_control%aux%VBias, &
                             SpinCL=smeagol_control%aux%SpinCL, &
                             NSlices=smeagol_control%aux%NSlices, &
                             TrCoeff=local_TrCoeff, &
                             CalcIETS=smeagol_control%aux%CalcIETS, &
                             ldos=local_ldos, &
                             idyn=0, & !
                             tmdskip=smeagol_control%aux%tmdskip, &
                             tmdsampling=smeagol_control%aux%tmdsampling)
      END IF

      ! Restore ovewriten EM.OrderN parameter
      IF (last) THEN
         smeagolglobal_negfon = negfon_saved
      END IF

      IF (PRESENT(rho_ao_kp)) THEN
         DO ispin = 1, nspin
            DO img = 1, nimages
               ! To be on a safe size, zeroize the new density matrix first. It is not actually needed.
               CALL dbcsr_set(rho_ao_kp(ispin, img)%matrix, 0.0_dp)
               matrix_kp_generic(img)%matrix => rho_ao_kp(ispin, img)%matrix
            END DO

            CALL convert_distributed_siesta_to_dbcsr(matrix_kp_generic, Dnew_csc_merged(:, ispin), &
                                                     siesta_struct_merged, para_env)
         END DO

         ! current-induced forces
         IF (smeagol_control%emforces) THEN
            Enew_csc_merged(:, :) = (1.0_dp/H_to_Ry)*Enew_csc_merged(:, :)

            DO ispin = 1, nspin
               DO img = 1, nimages
                  ! To be on a safe size, zeroize the new energy-density matrix first. It is not actually needed.
                  CALL dbcsr_set(matrix_w_kp(ispin, img)%matrix, 0.0_dp)
                  matrix_kp_generic(img)%matrix => matrix_w_kp(ispin, img)%matrix
               END DO

               CALL convert_distributed_siesta_to_dbcsr(matrix_kp_generic, Enew_csc_merged(:, ispin), &
                                                        siesta_struct_merged, para_env)
            END DO
         END IF
      END IF

      CALL destroy_communicators_negf()
      CALL emtrans_deallocate_global_arrays()

      DEALLOCATE (Dnew_csc_merged, Enew_csc_merged)
      DEALLOCATE (matrix_s_csc_merged, matrix_ks_csc_merged)

      CALL siesta_struct_release(siesta_struct_merged)

      IF (.NOT. do_kpoints) DEALLOCATE (cell_to_index)
#else
      CALL cp_abort(__LOCATION__, &
                    "CP2K was compiled with no SMEAGOL support.")
      MARK_USED(last)
      MARK_USED(iter)
      MARK_USED(rho_ao_kp)
      ! local variables
      MARK_USED(cell_to_index)
      MARK_USED(do_kpoints)
      MARK_USED(Dnew_csc_merged)
      MARK_USED(Enew_csc_merged)
      MARK_USED(img)
      MARK_USED(ispin)
      MARK_USED(local_ldos)
      MARK_USED(local_TrCoeff)
      MARK_USED(max_ij_cell_image)
      MARK_USED(matrix_kp_generic)
      MARK_USED(matrix_ks_csc_merged)
      MARK_USED(matrix_s_csc_merged)
      MARK_USED(md_step)
      MARK_USED(natoms)
      MARK_USED(negfon_saved)
      MARK_USED(nimages)
      MARK_USED(n_cell_images)
      MARK_USED(siesta_struct_merged)
      MARK_USED(structure_changed)
      MARK_USED(ucell)
#endif

      CALL timestop(handle)
   END SUBROUTINE run_smeagol_emtrans

END MODULE smeagol_interface
