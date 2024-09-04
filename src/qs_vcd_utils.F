!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_vcd_utils
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_result_methods,               ONLY: get_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE input_constants,                 ONLY: use_mom_ref_user
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_types,                  ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE particle_types,                  ONLY: particle_type
   USE qs_dcdr_utils,                   ONLY: dcdr_env_cleanup,&
                                              dcdr_env_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: vcd_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_local_moments_der_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_vcd_ao,                       ONLY: build_com_rpnl_r,&
                                              build_matrix_r_vhxc,&
                                              build_rcore_matrix,&
                                              build_rpnl_matrix,&
                                              build_tr_matrix
   USE string_utilities,                ONLY: xstring
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: vcd_env_cleanup, vcd_env_init
   PUBLIC :: vcd_read_restart, vcd_write_restart
   PUBLIC :: vcd_print

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vcd_utils'

   REAL(dp), DIMENSION(3, 3, 3), PARAMETER  :: Levi_Civita = RESHAPE((/ &
                                                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                          0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
                                             0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), (/3, 3, 3/))

CONTAINS

! *****************************************************************************
!> \brief Initialize the vcd environment
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE vcd_env_init(vcd_env, qs_env)
      TYPE(vcd_env_type), TARGET                         :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vcd_env_init'

      INTEGER                                            :: handle, i, idir, ispin, j, natom, &
                                                            nspins, output_unit, reference, &
                                                            unit_number
      LOGICAL                                            :: explicit
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, my_matrix_hr_1d
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: lr_section, vcd_section

      CALL timeset(routineN, handle)
      vcd_env%do_mfp = .FALSE.

      ! Set up the logger
      NULLIFY (logger, vcd_section, lr_section)
      logger => cp_get_default_logger()
      vcd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%VCD")
      vcd_env%output_unit = cp_print_key_unit_nr(logger, vcd_section, "PRINT%VCD", &
                                                 extension=".data", middle_name="vcd", log_filename=.FALSE., &
                                                 file_position="REWIND", file_status="REPLACE")

      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      output_unit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".linresLog")
      unit_number = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", extension=".linresLog")

      ! We can't run a NVPT/MFPT calculation without the coefficients dC/dR.
      CALL dcdr_env_init(vcd_env%dcdr_env, qs_env)
      ! vcd_env%dcdr_env%output_unit = vcd_env%output_unit

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(/,T20,A,/)") "*** Start NVPT/MFPT calculation ***"
      END IF

      ! Just to make sure. The memory requirements are tiny.
      CALL init_orbital_pointers(12)

      CALL section_vals_val_get(vcd_section, "DISTRIBUTED_ORIGIN", l_val=vcd_env%distributed_origin)
      CALL section_vals_val_get(vcd_section, "ORIGIN_DEPENDENT_MFP", l_val=vcd_env%origin_dependent_op_mfp)

      ! Reference point
      vcd_env%magnetic_origin = 0._dp
      vcd_env%spatial_origin = 0._dp
      ! Get the magnetic origin from the input
      CALL section_vals_val_get(vcd_section, "MAGNETIC_ORIGIN", i_val=reference)
      CALL section_vals_val_get(vcd_section, "MAGNETIC_ORIGIN_REFERENCE", explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(vcd_section, "MAGNETIC_ORIGIN_REFERENCE", r_vals=ref_point)
      ELSE
         IF (reference == use_mom_ref_user) &
            CPABORT("User-defined reference point should be given explicitly")
      END IF

      CALL get_reference_point(rpoint=vcd_env%magnetic_origin, qs_env=qs_env, &
                               reference=reference, &
                               ref_point=ref_point)

      ! Get the spatial origin from the input
      CALL section_vals_val_get(vcd_section, "SPATIAL_ORIGIN", i_val=reference)
      CALL section_vals_val_get(vcd_section, "SPATIAL_ORIGIN_REFERENCE", explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(vcd_section, "SPATIAL_ORIGIN_REFERENCE", r_vals=ref_point)
      ELSE
         IF (reference == use_mom_ref_user) &
            CPABORT("User-defined reference point should be given explicitly")
      END IF

      CALL get_reference_point(rpoint=vcd_env%spatial_origin, qs_env=qs_env, &
                               reference=reference, &
                               ref_point=ref_point)

      IF (vcd_env%distributed_origin .AND. ANY(vcd_env%magnetic_origin /= vcd_env%spatial_origin)) THEN
         CPWARN("The magnetic and spatial origins don't match")
         ! This is fine for NVP but will give unphysical results for MFP.
      END IF

      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,3F10.6)") &
         'The reference point is', vcd_env%dcdr_env%ref_point
      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,3F10.6)") &
         'The magnetic origin is', vcd_env%magnetic_origin
      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,3F10.6)") &
         'The velocity origin is', vcd_env%spatial_origin

      vcd_env%magnetic_origin_atom = vcd_env%magnetic_origin
      vcd_env%spatial_origin_atom = vcd_env%spatial_origin

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      sab_all=sab_all, &
                      sap_ppnl=sap_ppnl, &
                      particle_set=particle_set, &
                      matrix_ks=matrix_ks, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set)

      natom = SIZE(particle_set)
      nspins = dft_control%nspins

      ALLOCATE (vcd_env%apt_el_nvpt(3, 3, natom))
      ALLOCATE (vcd_env%apt_nuc_nvpt(3, 3, natom))
      ALLOCATE (vcd_env%apt_total_nvpt(3, 3, natom))
      ALLOCATE (vcd_env%aat_atom_nvpt(3, 3, natom))
      ALLOCATE (vcd_env%aat_atom_mfp(3, 3, natom))
      vcd_env%apt_el_nvpt = 0._dp
      vcd_env%apt_nuc_nvpt = 0._dp
      vcd_env%apt_total_nvpt = 0._dp
      vcd_env%aat_atom_nvpt = 0._dp
      vcd_env%aat_atom_mfp = 0._dp

      ALLOCATE (vcd_env%dCV(nspins))
      ALLOCATE (vcd_env%dCV_prime(nspins))
      ALLOCATE (vcd_env%op_dV(nspins))
      ALLOCATE (vcd_env%op_dB(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_create(vcd_env%dCV(ispin), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
         CALL cp_fm_create(vcd_env%dCV_prime(ispin), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
         CALL cp_fm_create(vcd_env%op_dV(ispin), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
         CALL cp_fm_create(vcd_env%op_dB(ispin), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
      END DO

      ALLOCATE (vcd_env%dCB(3))
      ALLOCATE (vcd_env%dCB_prime(3))
      DO i = 1, 3
         CALL cp_fm_create(vcd_env%dCB(i), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
         CALL cp_fm_create(vcd_env%dCB_prime(i), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
      END DO

      ! DBCSR matrices
      CALL dbcsr_allocate_matrix_set(vcd_env%moments_der, 9, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%moments_der_right, 9, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%moments_der_left, 9, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_difdip2, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_dSdV, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_dSdB, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_hxc_dsdv, nspins)

      CALL dbcsr_allocate_matrix_set(vcd_env%hcom, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_rcomr, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_rrcom, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_dcom, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_hr, nspins, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_rh, nspins, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_drpnl, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%dipvel_ao, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%dipvel_ao_delta, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_rxrv, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_r_rxvr, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_rxvr_r, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_r_doublecom, 3, 3)

      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_nosym_temp_33, 3, 3)
      CALL dbcsr_allocate_matrix_set(vcd_env%matrix_nosym_temp2_33, 3, 3)
      DO i = 1, 9 ! x, y, z, xx, xy, xz, yy, yz, zz
         DO idir = 1, 3 ! d/dx, d/dy, d/dz
            CALL dbcsr_init_p(vcd_env%moments_der(i, idir)%matrix)
            CALL dbcsr_init_p(vcd_env%moments_der_right(i, idir)%matrix)
            CALL dbcsr_init_p(vcd_env%moments_der_left(i, idir)%matrix)

            CALL dbcsr_create(vcd_env%moments_der(i, idir)%matrix, template=matrix_ks(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(vcd_env%moments_der(i, idir)%matrix, sab_orb)
            CALL dbcsr_set(vcd_env%moments_der(i, idir)%matrix, 0.0_dp)

            ! And the ones which will be multiplied by delta_(mu/nu)
            CALL dbcsr_copy(vcd_env%moments_der_right(i, idir)%matrix, &
                            vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_copy(vcd_env%moments_der_left(i, idir)%matrix, &
                            vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
         END DO
      END DO

      DO i = 1, 3
         DO j = 1, 3
            CALL dbcsr_init_p(vcd_env%matrix_difdip2(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_difdip2(i, j)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_set(vcd_env%matrix_difdip2(i, j)%matrix, 0.0_dp)

            CALL dbcsr_init_p(vcd_env%matrix_nosym_temp_33(i, j)%matrix)
            CALL dbcsr_create(vcd_env%matrix_nosym_temp_33(i, j)%matrix, template=matrix_ks(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(vcd_env%matrix_nosym_temp_33(i, j)%matrix, sab_all)
            CALL dbcsr_set(vcd_env%matrix_nosym_temp_33(i, j)%matrix, 0._dp)

            CALL dbcsr_init_p(vcd_env%matrix_nosym_temp2_33(i, j)%matrix)
            CALL dbcsr_create(vcd_env%matrix_nosym_temp2_33(i, j)%matrix, template=matrix_ks(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(vcd_env%matrix_nosym_temp2_33(i, j)%matrix, sab_all)
            CALL dbcsr_set(vcd_env%matrix_nosym_temp2_33(i, j)%matrix, 0._dp)

         END DO
         CALL dbcsr_init_p(vcd_env%matrix_dSdV(i)%matrix)
         CALL dbcsr_copy(vcd_env%matrix_dSdV(i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
         CALL dbcsr_init_p(vcd_env%matrix_dSdB(i)%matrix)
         CALL dbcsr_copy(vcd_env%matrix_dSdB(i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
      END DO

      DO ispin = 1, nspins
         CALL dbcsr_init_p(vcd_env%matrix_hxc_dsdv(ispin)%matrix)
         CALL dbcsr_copy(vcd_env%matrix_hxc_dsdv(ispin)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
      END DO

      ! Things for op_dV
      ! lin_mom
      DO i = 1, 3
         CALL dbcsr_init_p(vcd_env%dipvel_ao(i)%matrix)
         CALL dbcsr_copy(vcd_env%dipvel_ao(i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)

         CALL dbcsr_init_p(vcd_env%dipvel_ao_delta(i)%matrix)
         CALL dbcsr_copy(vcd_env%dipvel_ao_delta(i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
      END DO

      ! [V, r]
      DO i = 1, 3
         CALL dbcsr_init_p(vcd_env%hcom(i)%matrix)
         CALL dbcsr_create(vcd_env%hcom(i)%matrix, template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(vcd_env%hcom(i)%matrix, sab_orb)

         CALL dbcsr_init_p(vcd_env%matrix_rxrv(i)%matrix)
         CALL dbcsr_create(vcd_env%matrix_rxrv(i)%matrix, template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(vcd_env%matrix_rxrv(i)%matrix, sab_orb)

         DO j = 1, 3
            CALL dbcsr_init_p(vcd_env%matrix_rcomr(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_rcomr(i, j)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_init_p(vcd_env%matrix_rrcom(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_rrcom(i, j)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_init_p(vcd_env%matrix_dcom(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_dcom(i, j)%matrix, matrix_ks(1)%matrix)
            CALL dbcsr_set(vcd_env%matrix_dcom(i, j)%matrix, 0._dp)

            CALL dbcsr_init_p(vcd_env%matrix_r_rxvr(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_r_rxvr(i, j)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_set(vcd_env%matrix_r_rxvr(i, j)%matrix, 0._dp)

            CALL dbcsr_init_p(vcd_env%matrix_rxvr_r(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_rxvr_r(i, j)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_set(vcd_env%matrix_rxvr_r(i, j)%matrix, 0._dp)

            CALL dbcsr_init_p(vcd_env%matrix_r_doublecom(i, j)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_r_doublecom(i, j)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_set(vcd_env%matrix_r_doublecom(i, j)%matrix, 0._dp)
         END DO
      END DO

      ! matrix_hr: nonsymmetric dbcsr matrix
      DO ispin = 1, nspins
         DO i = 1, 3
            CALL dbcsr_init_p(vcd_env%matrix_hr(ispin, i)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_hr(ispin, i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)

            CALL dbcsr_init_p(vcd_env%matrix_rh(ispin, i)%matrix)
            CALL dbcsr_copy(vcd_env%matrix_rh(ispin, i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
         END DO
      END DO

      ! drpnl for the operator
      DO i = 1, 3
         CALL dbcsr_init_p(vcd_env%matrix_drpnl(i)%matrix)
         CALL dbcsr_copy(vcd_env%matrix_drpnl(i)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
      END DO

      ! NVP matrices
      ! hr matrices
      my_matrix_hr_1d => vcd_env%matrix_hr(1, 1:3)
      CALL build_rpnl_matrix(my_matrix_hr_1d, qs_kind_set, particle_set, sab_all, sap_ppnl, &
                             dft_control%qs_control%eps_ppnl, cell, [0._dp, 0._dp, 0._dp], &
                             direction_Or=.TRUE.)
      CALL build_tr_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, &
                           direction_Or=.TRUE., rc=[0._dp, 0._dp, 0._dp])
      CALL build_rcore_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, [0._dp, 0._dp, 0._dp])
      CALL build_matrix_r_vhxc(vcd_env%matrix_hr, qs_env, [0._dp, 0._dp, 0._dp])

      my_matrix_hr_1d => vcd_env%matrix_rh(1, 1:3)
      CALL build_rpnl_matrix(my_matrix_hr_1d, qs_kind_set, particle_set, sab_all, sap_ppnl, &
                             dft_control%qs_control%eps_ppnl, cell, [0._dp, 0._dp, 0._dp], &
                             direction_Or=.FALSE.)
      CALL build_tr_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, &
                           direction_Or=.FALSE., rc=[0._dp, 0._dp, 0._dp])
      CALL build_rcore_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, [0._dp, 0._dp, 0._dp])
      CALL build_matrix_r_vhxc(vcd_env%matrix_rh, qs_env, [0._dp, 0._dp, 0._dp])

      ! commutator terms
      ! - [V, r]
      CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                            particle_set, cell=cell, matrix_rv=vcd_env%hcom)
      ! <[V, r] * r> and <r * [V, r]>
      CALL build_com_rpnl_r(vcd_env%matrix_rcomr, qs_kind_set, sab_all, sap_ppnl, &
                            dft_control%qs_control%eps_ppnl, particle_set, cell, .TRUE.)
      CALL build_com_rpnl_r(vcd_env%matrix_rrcom, qs_kind_set, sab_all, sap_ppnl, &
                            dft_control%qs_control%eps_ppnl, particle_set, cell, .FALSE.)

      ! lin_mom
      CALL build_lin_mom_matrix(qs_env, vcd_env%dipvel_ao)

      ! AAT
      ! The moments are set to zero and then recomputed in the routine.
      CALL build_local_moments_der_matrix(qs_env, moments_der=vcd_env%moments_der, &
                                          nmoments_der=2, nmoments=0, ref_point=[0._dp, 0._dp, 0._dp])

      ! PP terms
      CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                            particle_set, matrix_rxrv=vcd_env%matrix_rxrv, ref_point=[0._dp, 0._dp, 0._dp], &
                            cell=cell)

      CALL build_com_mom_nl(qs_kind_set, sab_all, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                            particle_set, ref_point=[0._dp, 0._dp, 0._dp], cell=cell, &
                            matrix_r_rxvr=vcd_env%matrix_r_rxvr)

      CALL build_com_mom_nl(qs_kind_set, sab_all, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                            particle_set, ref_point=[0._dp, 0._dp, 0._dp], cell=cell, &
                            matrix_rxvr_r=vcd_env%matrix_rxvr_r)

      ! Done with NVP matrices

      CALL cp_print_key_finished_output(output_unit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE vcd_env_init

! *****************************************************************************
!> \brief Deallocate the vcd environment
!> \param qs_env ...
!> \param vcd_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE vcd_env_cleanup(qs_env, vcd_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(vcd_env_type)                                 :: vcd_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vcd_env_cleanup'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! We can't run a NVPT/MFPT calculation without the coefficients dC/dR.
      CALL dcdr_env_cleanup(qs_env, vcd_env%dcdr_env)

      DEALLOCATE (vcd_env%apt_el_nvpt)
      DEALLOCATE (vcd_env%apt_nuc_nvpt)
      DEALLOCATE (vcd_env%apt_total_nvpt)
      DEALLOCATE (vcd_env%aat_atom_nvpt)
      DEALLOCATE (vcd_env%aat_atom_mfp)

      CALL cp_fm_release(vcd_env%dCV)
      CALL cp_fm_release(vcd_env%dCV_prime)
      CALL cp_fm_release(vcd_env%op_dV)
      CALL cp_fm_release(vcd_env%op_dB)

      CALL cp_fm_release(vcd_env%dCB)
      CALL cp_fm_release(vcd_env%dCB_prime)

      ! DBCSR matrices
      ! Probably, the memory requirements could be reduced by quite a bit
      ! by not storing each term in its own set of matrices.
      ! On the other hand, the memory bottleneck is usually the numerical
      ! integration grid.
      CALL dbcsr_deallocate_matrix_set(vcd_env%moments_der)
      CALL dbcsr_deallocate_matrix_set(vcd_env%moments_der_right)
      CALL dbcsr_deallocate_matrix_set(vcd_env%moments_der_left)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_difdip2)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_dSdV)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_dSdB)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_hxc_dsdv)
      CALL dbcsr_deallocate_matrix_set(vcd_env%hcom)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_rcomr)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_rrcom)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_dcom)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_hr)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_rh)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_drpnl)
      CALL dbcsr_deallocate_matrix_set(vcd_env%dipvel_ao)
      CALL dbcsr_deallocate_matrix_set(vcd_env%dipvel_ao_delta)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_rxrv)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_r_rxvr)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_rxvr_r)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_r_doublecom)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_nosym_temp_33)
      CALL dbcsr_deallocate_matrix_set(vcd_env%matrix_nosym_temp2_33)
      CALL timestop(handle)

   END SUBROUTINE vcd_env_cleanup

! **************************************************************************************************
!> \brief Copied from linres_read_restart
!> \param qs_env ...
!> \param linres_section ...
!> \param vec ...
!> \param lambda ...
!> \param beta ...
!> \param tag ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE vcd_read_restart(qs_env, linres_section, vec, lambda, beta, tag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: linres_section
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: vec
      INTEGER, INTENT(IN)                                :: lambda, beta
      CHARACTER(LEN=*)                                   :: tag

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vcd_read_restart'

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: my_middle
      INTEGER :: beta_tmp, handle, i, i_block, ia, ie, iostat, iounit, ispin, j, lambda_tmp, &
         max_block, n_rep_val, nao, nao_tmp, nmo, nmo_tmp, nspins, nspins_tmp, rst_unit
      LOGICAL                                            :: file_exists
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: print_key

      file_exists = .FALSE.

      CALL timeset(routineN, handle)

      NULLIFY (mos, para_env, logger, print_key, vecbuffer)
      logger => cp_get_default_logger()

      iounit = cp_print_key_unit_nr(logger, linres_section, &
                                    "PRINT%PROGRAM_RUN_INFO", extension=".Log")

      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env, &
                      mos=mos)

      nspins = SIZE(mos)

      rst_unit = -1
      IF (para_env%is_source()) THEN
         CALL section_vals_val_get(linres_section, "WFN_RESTART_FILE_NAME", &
                                   n_rep_val=n_rep_val)

         CALL XSTRING(tag, ia, ie)
         my_middle = "RESTART-"//tag(ia:ie)//TRIM("-")//TRIM(ADJUSTL(cp_to_string(beta))) &
                     //TRIM("-")//TRIM(ADJUSTL(cp_to_string(lambda)))

         IF (n_rep_val > 0) THEN
            CALL section_vals_val_get(linres_section, "WFN_RESTART_FILE_NAME", c_val=filename)
            CALL xstring(filename, ia, ie)
            filename = filename(ia:ie)//TRIM(my_middle)//".lr"
         ELSE
            ! try to read from the filename that is generated automatically from the printkey
            print_key => section_vals_get_subs_vals(linres_section, "PRINT%RESTART")
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      extension=".lr", middle_name=TRIM(my_middle), my_local=.FALSE.)
         END IF
         INQUIRE (FILE=filename, exist=file_exists)
         !
         ! open file
         IF (file_exists) THEN
            CALL open_file(file_name=TRIM(filename), &
                           file_action="READ", &
                           file_form="UNFORMATTED", &
                           file_position="REWIND", &
                           file_status="OLD", &
                           unit_number=rst_unit)

            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Reading response wavefunctions from the restart file <"//TRIM(ADJUSTL(filename))//">"
         ELSE
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file  <"//TRIM(ADJUSTL(filename))//"> not found"
         END IF
      END IF

      CALL para_env%bcast(file_exists)

      IF (file_exists) THEN

         CALL get_mo_set(mos(1), mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_block=max_block)

         ALLOCATE (vecbuffer(nao, max_block))
         !
         ! read headers
         IF (rst_unit > 0) READ (rst_unit, IOSTAT=iostat) lambda_tmp, beta_tmp, nspins_tmp, nao_tmp
         CALL para_env%bcast(iostat)

         CALL para_env%bcast(beta_tmp)
         CALL para_env%bcast(lambda_tmp)
         CALL para_env%bcast(nspins_tmp)
         CALL para_env%bcast(nao_tmp)

         ! check that the number nao, nmo and nspins are
         ! the same as in the current mos
         IF (nspins_tmp .NE. nspins) THEN
            CPABORT("nspins not consistent")
         END IF
         IF (nao_tmp .NE. nao) CPABORT("nao not consistent")
         ! check that it's the right file
         ! the same as in the current mos
         IF (lambda_tmp .NE. lambda) CPABORT("lambda not consistent")
         IF (beta_tmp .NE. beta) CPABORT("beta not consistent")
         !
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, ncol_global=nmo)
            !
            IF (rst_unit > 0) READ (rst_unit) nmo_tmp
            CALL para_env%bcast(nmo_tmp)
            IF (nmo_tmp .NE. nmo) CPABORT("nmo not consistent")
            !
            ! read the response
            DO i = 1, nmo, MAX(max_block, 1)
               i_block = MIN(max_block, nmo - i + 1)
               DO j = 1, i_block
                  IF (rst_unit > 0) READ (rst_unit) vecbuffer(1:nao, j)
               END DO
               CALL para_env%bcast(vecbuffer)
               CALL cp_fm_set_submatrix(vec(ispin), vecbuffer, 1, i, nao, i_block)
            END DO
         END DO

         IF (iostat /= 0) THEN
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file <"//TRIM(ADJUSTL(filename))//"> not found"
         END IF

         DEALLOCATE (vecbuffer)

      END IF

      IF (para_env%is_source()) THEN
         IF (file_exists) CALL close_file(unit_number=rst_unit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE vcd_read_restart

! **************************************************************************************************
!> \brief Copied from linres_write_restart
!> \param qs_env ...
!> \param linres_section ...
!> \param vec ...
!> \param lambda ...
!> \param beta ...
!> \param tag ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE vcd_write_restart(qs_env, linres_section, vec, lambda, beta, tag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: linres_section
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: vec
      INTEGER, INTENT(IN)                                :: lambda, beta
      CHARACTER(LEN=*)                                   :: tag

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vcd_write_restart'

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: my_middle, my_pos, my_status
      INTEGER                                            :: handle, i, i_block, ia, ie, iounit, &
                                                            ispin, j, max_block, nao, nmo, nspins, &
                                                            rst_unit
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, mo_coeff, mos, para_env, print_key, vecbuffer)

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, linres_section, "PRINT%RESTART", &
                                           used_print_key=print_key), &
                cp_p_file)) THEN

         iounit = cp_print_key_unit_nr(logger, linres_section, &
                                       "PRINT%PROGRAM_RUN_INFO", extension=".Log")

         CALL get_qs_env(qs_env=qs_env, &
                         mos=mos, &
                         para_env=para_env)

         nspins = SIZE(mos)

         my_status = "REPLACE"
         my_pos = "REWIND"
         CALL XSTRING(tag, ia, ie)
         my_middle = "RESTART-"//tag(ia:ie)//TRIM("-")//TRIM(ADJUSTL(cp_to_string(beta))) &
                     //TRIM("-")//TRIM(ADJUSTL(cp_to_string(lambda)))
         rst_unit = cp_print_key_unit_nr(logger, linres_section, "PRINT%RESTART", &
                                         extension=".lr", middle_name=TRIM(my_middle), file_status=TRIM(my_status), &
                                         file_position=TRIM(my_pos), file_action="WRITE", file_form="UNFORMATTED")

         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".lr", middle_name=TRIM(my_middle), my_local=.FALSE.)

         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "LINRES| Writing response functions to the restart file <"//TRIM(ADJUSTL(filename))//">"
         END IF

         !
         ! write data to file
         ! use the scalapack block size as a default for buffering columns
         CALL get_mo_set(mos(1), mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_block=max_block)
         ALLOCATE (vecbuffer(nao, max_block))

         IF (rst_unit > 0) WRITE (rst_unit) lambda, beta, nspins, nao

         DO ispin = 1, nspins
            CALL cp_fm_get_info(vec(ispin), ncol_global=nmo)

            IF (rst_unit > 0) WRITE (rst_unit) nmo

            DO i = 1, nmo, MAX(max_block, 1)
               i_block = MIN(max_block, nmo - i + 1)
               CALL cp_fm_get_submatrix(vec(ispin), vecbuffer, 1, i, nao, i_block)
               ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
               ! to old ones, and in cases where max_block is different between runs, as might happen during
               ! restarts with a different number of CPUs
               DO j = 1, i_block
                  IF (rst_unit > 0) WRITE (rst_unit) vecbuffer(1:nao, j)
               END DO
            END DO
         END DO

         DEALLOCATE (vecbuffer)

         CALL cp_print_key_finished_output(rst_unit, logger, linres_section, &
                                           "PRINT%RESTART")
      END IF

      CALL timestop(handle)

   END SUBROUTINE vcd_write_restart

! **************************************************************************************************
!> \brief Print the APTs, AATs, and sum rules
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE vcd_print(vcd_env, qs_env)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'vcd_print'

      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: alpha, beta, delta, gamma, handle, i, l, &
                                                            lambda, natom, nsubset, output_unit
      REAL(dp)                                           :: mean, standard_deviation, &
                                                            standard_deviation_sum
      REAL(dp), DIMENSION(:, :, :), POINTER              :: apt_el_dcdr, apt_el_nvpt, apt_nuc_dcdr, &
                                                            apt_nuc_nvpt, apt_total_dcdr, &
                                                            apt_total_nvpt
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: apt_center_dcdr, apt_subset_dcdr
      REAL(kind=dp), DIMENSION(3, 3)                     :: sum_rule_0, sum_rule_0_second, &
                                                            sum_rule_1, sum_rule_2, &
                                                            sum_rule_2_second, sum_rule_3_mfp, &
                                                            sum_rule_3_second
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: vcd_section

      CALL timeset(routineN, handle)

      NULLIFY (logger)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      vcd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%VCD")

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, molecule_set=molecule_set)
      natom = SIZE(particle_set)
      nsubset = SIZE(molecule_set)

      apt_el_dcdr => vcd_env%dcdr_env%apt_el_dcdr
      apt_nuc_dcdr => vcd_env%dcdr_env%apt_nuc_dcdr
      apt_total_dcdr => vcd_env%dcdr_env%apt_total_dcdr
      apt_subset_dcdr => vcd_env%dcdr_env%apt_el_dcdr_per_subset
      apt_center_dcdr => vcd_env%dcdr_env%apt_el_dcdr_per_center

      apt_el_nvpt => vcd_env%apt_el_nvpt
      apt_nuc_nvpt => vcd_env%apt_nuc_nvpt
      apt_total_nvpt => vcd_env%apt_total_nvpt

      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A)") &
         'APT | Write the final APT matrix per atom (Position perturbation)'
      DO l = 1, natom
         IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,I3,A,F15.6)") &
            'APT | Atom', l, ' - GAPT ', &
            (apt_total_dcdr(1, 1, l) &
             + apt_total_dcdr(2, 2, l) &
             + apt_total_dcdr(3, 3, l))/3._dp
         DO i = 1, 3
            IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,F15.6,F15.6,F15.6)") "APT | ", apt_total_dcdr(i, :, l)
         END DO
      END DO

      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A)") &
         'NVP | Write the final APT matrix per atom (Velocity perturbation)'
      DO l = 1, natom
         IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,I3,A,F15.6)") &
            'NVP | Atom', l, ' - GAPT ', &
            (apt_total_nvpt(1, 1, l) &
             + apt_total_nvpt(2, 2, l) &
             + apt_total_nvpt(3, 3, l))/3._dp
         DO i = 1, 3
            IF (vcd_env%output_unit > 0) &
               WRITE (vcd_env%output_unit, "(A,F15.6,F15.6,F15.6)") &
               "NVP | ", apt_total_nvpt(i, :, l)
         END DO
      END DO

      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A)") &
         'NVP | Write the final AAT matrix per atom (Velocity perturbation)'
      DO l = 1, natom
         IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,I3)") &
            'NVP | Atom', l
         DO i = 1, 3
            IF (vcd_env%output_unit > 0) &
               WRITE (vcd_env%output_unit, "(A,F15.6,F15.6,F15.6)") &
               "NVP | ", vcd_env%aat_atom_nvpt(i, :, l)
         END DO
      END DO

      IF (vcd_env%do_mfp) THEN
         IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A)") &
            'MFP | Write the final AAT matrix per atom (Magnetic Field perturbation)'
         DO l = 1, natom
            IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,I3)") &
               'MFP | Atom', l
            DO i = 1, 3
               IF (vcd_env%output_unit > 0) &
                  WRITE (vcd_env%output_unit, "(A,F15.6,F15.6,F15.6)") &
                  "MFP | ", vcd_env%aat_atom_mfp(i, :, l)
            END DO
         END DO
      END IF

      ! Get the dipole
      CALL get_qs_env(qs_env, results=results)
      description = "[DIPOLE]"
      CALL get_results(results=results, description=description, values=vcd_env%dcdr_env%dipole_pos(1:3))

      ! Sum rules [for all alpha, beta]
      sum_rule_0 = 0._dp
      sum_rule_1 = 0._dp
      sum_rule_2 = 0._dp
      sum_rule_0_second = 0._dp
      sum_rule_2_second = 0._dp
      sum_rule_3_second = 0._dp
      sum_rule_3_mfp = 0._dp
      standard_deviation = 0._dp
      standard_deviation_sum = 0._dp

      DO alpha = 1, 3
         DO beta = 1, 3
            ! 0: sum_lambda apt(alpha, beta, lambda)
            DO lambda = 1, natom
               sum_rule_0(alpha, beta) = sum_rule_0(alpha, beta) &
                                         + apt_total_dcdr(alpha, beta, lambda)
               sum_rule_0_second(alpha, beta) = sum_rule_0_second(alpha, beta) &
                                                + apt_total_nvpt(alpha, beta, lambda)
            END DO

            ! 1: sum_gamma epsilon_(alpha beta gamma) mu_gamma
            DO gamma = 1, 3
               sum_rule_1(alpha, beta) = sum_rule_1(alpha, beta) &
                                         + Levi_Civita(alpha, beta, gamma)*vcd_env%dcdr_env%dipole_pos(gamma)
            END DO

            ! 2: sum_(lambda gamma delta) R^lambda_gamma apt(delta, alpha, lambda)
            DO lambda = 1, natom
               DO gamma = 1, 3
                  DO delta = 1, 3
                     sum_rule_2(alpha, beta) = sum_rule_2(alpha, beta) &
                                               + Levi_Civita(beta, gamma, delta) &
                                               *particle_set(lambda)%r(gamma) &
                                               *apt_total_dcdr(delta, alpha, lambda)
                     sum_rule_2_second(alpha, beta) = sum_rule_2_second(alpha, beta) &
                                                      + Levi_Civita(beta, gamma, delta) &
                                                      *particle_set(lambda)%r(gamma) &
                                                      *apt_total_nvpt(delta, alpha, lambda)
                  END DO
               END DO
            END DO

            ! 3: 2c * sum_lambda aat(alpha, beta, lambda)
            DO lambda = 1, natom
               sum_rule_3_second(alpha, beta) = sum_rule_3_second(alpha, beta) &
                                                + vcd_env%aat_atom_nvpt(alpha, beta, lambda)
               ! + 2._dp*c_light_au*vcd_env%aat_atom_nvpt(alpha, beta, lambda)
            END DO

            IF (vcd_env%do_mfp) THEN
               ! 3: 2c * sum_lambda aat(alpha, beta, lambda)
               DO lambda = 1, natom
                  sum_rule_3_mfp(alpha, beta) = sum_rule_3_mfp(alpha, beta) &
                                                + vcd_env%aat_atom_mfp(alpha, beta, lambda)
               END DO
            END IF

         END DO ! beta
      END DO   ! alpha

      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A)") "APT | Position perturbation sum rules"
      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(A,T19,A,T35,A,T50,A,T65,A)") &
         "APT |", " Total APT", "Dipole", "R * APT", "AAT"
      standard_deviation_sum = 0._dp
      DO alpha = 1, 3
         DO beta = 1, 3
            mean = (sum_rule_1(alpha, beta) + sum_rule_2(alpha, beta) + sum_rule_3_mfp(alpha, beta))/3
            standard_deviation = &
               SQRT((sum_rule_1(alpha, beta)**2 + sum_rule_2(alpha, beta)**2 + sum_rule_3_mfp(alpha, beta)**2)/3 &
                    - mean**2)
            standard_deviation_sum = standard_deviation_sum + standard_deviation

            IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, &
                                                "(A,I3,I3,F15.6,F15.6,F15.6,F15.6,F15.6)") &
               "APT | ", &
               alpha, beta, &
               sum_rule_0(alpha, beta), &
               sum_rule_1(alpha, beta), &
               sum_rule_2(alpha, beta), &
               sum_rule_3_mfp(alpha, beta), &
               standard_deviation
         END DO
      END DO
      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(T73,F15.6)") standard_deviation_sum

      IF (vcd_env%output_unit > 0) THEN
         WRITE (vcd_env%output_unit, "(A)") "NVP | Velocity perturbation sum rules"
         WRITE (vcd_env%output_unit, "(A,T19,A,T35,A,T50,A,T65,A)") "NVP |", " Total APT", "Dipole", "R * APT", "AAT"
      END IF

      standard_deviation_sum = 0._dp
      DO alpha = 1, 3
         DO beta = 1, 3
            mean = (sum_rule_1(alpha, beta) + sum_rule_2_second(alpha, beta) + sum_rule_3_second(alpha, beta))/3
            standard_deviation = &
               SQRT((sum_rule_1(alpha, beta)**2 + sum_rule_2_second(alpha, beta)**2 + sum_rule_3_second(alpha, beta)**2)/3 &
                    - mean**2)
            standard_deviation_sum = standard_deviation_sum + standard_deviation
            IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, &
                                                "(A,I3,I3,F15.6,F15.6,F15.6,F15.6,F15.6)") &
               "NVP | ", &
               alpha, &
               beta, &
               sum_rule_0_second(alpha, beta), &
               sum_rule_1(alpha, beta), &
               sum_rule_2_second(alpha, beta), &
               sum_rule_3_second(alpha, beta), &
               standard_deviation
         END DO
      END DO
      IF (vcd_env%output_unit > 0) WRITE (vcd_env%output_unit, "(T73,F15.6)") standard_deviation_sum

      CALL timestop(handle)
   END SUBROUTINE vcd_print

END MODULE qs_vcd_utils
