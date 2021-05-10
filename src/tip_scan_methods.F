!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE tip_scan_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: silent_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_cube_to_pw
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_bcast
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release,&
                                              pw_grid_setup
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_structure_factor,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              mo_set_p_type,&
                                              reassign_allocated_mos
   USE qs_scf,                          ONLY: scf
   USE tip_scan_types,                  ONLY: read_scanning_section,&
                                              release_scanning_type,&
                                              scanning_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tip_scan_methods'

   PUBLIC :: tip_scanning

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Perform tip scanning calculation.
!> \param qs_env  Quickstep environment
!>        input_section  Tip Scan Section
!> \param input_section ...
!> \par History
!>    * 05.2021 created [JGH]
! **************************************************************************************************
   SUBROUTINE tip_scanning(qs_env, input_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tip_scanning'

      CHARACTER(LEN=default_string_length)               :: cname
      INTEGER                                            :: handle, iounit, iscan, iset, nscan, &
                                                            nset, plevel, tsteps
      LOGICAL                                            :: do_tip_scan, expot, scf_converged
      REAL(KIND=dp), DIMENSION(3)                        :: rpos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mos_ref
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), POINTER                           :: sf, vee, vref, vtip
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(scanning_type)                                :: scan_env

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      CALL section_vals_val_get(input_section, "_SECTION_PARAMETERS_", l_val=do_tip_scan)
      IF (do_tip_scan) THEN
         iounit = cp_logger_get_default_io_unit(logger)
         cname = logger%iter_info%project_name
         logger%iter_info%project_name = logger%iter_info%project_name//"+TIP_SCAN"
         plevel = logger%iter_info%print_level
         logger%iter_info%print_level = silent_print_level

         IF (iounit > 0) THEN
            WRITE (iounit, "(T2,A)") "TIP SCAN| Perform a Tip Scanning Calculation"
         END IF

         ! read the input section
         CALL read_scanning_section(scan_env, input_section)
         ! read tip potential file
         CALL read_tip_file(qs_env, scan_env)

         CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env, &
                         dft_control=dft_control)
         expot = dft_control%apply_external_potential
         dft_control%apply_external_potential = .TRUE.
         NULLIFY (vee, vtip)
         IF (expot) THEN
            ! save external potential
            CALL get_qs_env(qs_env, vee=vee)
         END IF

         ! scratch memory for tip potentials and structure factor
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         ALLOCATE (vtip)
         CALL pw_pool_create_pw(auxbas_pw_pool, vtip%pw, use_data=REALDATA3D)
         vtip%pw%in_space = REALSPACE
         CALL pw_zero(vtip%pw)
         NULLIFY (sf, vref)
         ALLOCATE (sf, vref)
         CALL pw_pool_create_pw(auxbas_pw_pool, vref%pw, use_data=COMPLEXDATA1D)
         vref%pw%in_space = RECIPROCALSPACE
         CALL pw_zero(vref%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, sf%pw, use_data=COMPLEXDATA1D)
         sf%pw%in_space = RECIPROCALSPACE

         ! get the reference tip potential and store it in reciprocal space (vref)
         CALL pw_transfer(scan_env%tip_pw_g, vref%pw)

         ! store reference MOs
         CALL get_qs_env(qs_env, mos=mos)
         nset = SIZE(mos)
         ALLOCATE (mos_ref(nset))
         DO iset = 1, nset
            CALL duplicate_mo_set(mos_ref(iset)%mo_set, mos(iset)%mo_set)
         END DO

         nscan = scan_env%num_scan_points
         IF (iounit > 0) THEN
            WRITE (iounit, "(T2,A,T74,I7)") "TIP SCAN| Number of scanning points ", nscan
            WRITE (iounit, "(T2,A)") "TIP SCAN| Start scanning ..."
         END IF

         DO iscan = 1, nscan
            IF (iounit > 0) THEN
               WRITE (iounit, "(T2,A,I7)", advance="NO") "TIP SCAN| Scan point ", iscan
            END IF

            ! shift the reference tip potential
            rpos(1:3) = scan_env%tip_pos(1:3, iscan) - scan_env%ref_point(1:3)
            CALL shift_tip_potential(vref, sf, vtip, rpos)
            ! set the external potential
            IF (ASSOCIATED(vee)) THEN
               CALL pw_axpy(vee%pw, vtip%pw, alpha=1.0_dp)
            END IF
            CALL set_ks_env(ks_env, vee=vtip)

            ! reset MOs
            CALL get_qs_env(qs_env, mos=mos)
            DO iset = 1, nset
               CALL reassign_allocated_mos(mos(iset)%mo_set, mos_ref(iset)%mo_set)
            END DO

            ! Calculate electronic structure
            CALL scf(qs_env, has_converged=scf_converged, total_scf_steps=tsteps)

            IF (iounit > 0) THEN
               IF (scf_converged) THEN
                  WRITE (iounit, "(T25,A,I4,A)") "SCF converged in ", tsteps, " steps"
               ELSE
                  WRITE (iounit, "(T31,A)") "SCF did not converge!"
               END IF
            END IF
         END DO
         CALL release_scanning_type(scan_env)

         IF (iounit > 0) THEN
            WRITE (iounit, "(T2,A)") "TIP SCAN| ... end scanning"
         END IF
         dft_control%apply_external_potential = expot
         IF (expot) THEN
            ! restore vee
            CALL set_ks_env(ks_env, vee=vee)
         ELSE
            NULLIFY (vee)
            CALL set_ks_env(ks_env, vee=vee)
         END IF
         CALL pw_pool_give_back_pw(auxbas_pw_pool, vtip%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, vref%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, sf%pw)
         DEALLOCATE (vtip, sf, vref)

         logger%iter_info%print_level = plevel
         logger%iter_info%project_name = cname

         ! reset MOs
         CALL get_qs_env(qs_env, mos=mos)
         DO iset = 1, nset
            CALL reassign_allocated_mos(mos(iset)%mo_set, mos_ref(iset)%mo_set)
            CALL deallocate_mo_set(mos_ref(iset)%mo_set)
         END DO
         DEALLOCATE (mos_ref)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tip_scanning

! **************************************************************************************************
!> \brief Shift tip potential in reciprocal space
!> \param vref ...
!> \param sf ...
!> \param vtip ...
!> \param rpos ...
!> \par History
!>    * 05.2021 created [JGH]
! **************************************************************************************************
   SUBROUTINE shift_tip_potential(vref, sf, vtip, rpos)

      TYPE(pw_p_type), POINTER                           :: vref, sf, vtip
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rpos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'shift_tip_potential'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL pw_structure_factor(sf%pw, rpos)
      sf%pw%cc = vref%pw%cc*sf%pw%cc
      CALL pw_transfer(sf%pw, vtip%pw)

      CALL timestop(handle)

   END SUBROUTINE shift_tip_potential

! **************************************************************************************************
!> \brief Read tip potential from cube file. Allow any spacing and cell size
!> \param qs_env ...
!> \param scan_env ...
!> \par History
!>    * 05.2021 created [JGH]
! **************************************************************************************************
   SUBROUTINE read_tip_file(qs_env, scan_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(scanning_type), INTENT(INOUT)                 :: scan_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_tip_file'

      INTEGER                                            :: extunit, i, nat
      INTEGER, DIMENSION(3)                              :: npts
      REAL(KIND=dp)                                      :: scaling
      REAL(KIND=dp), DIMENSION(3)                        :: rdum
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dcell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_type), POINTER                             :: pw

      CALL get_qs_env(qs_env, para_env=para_env)

      IF (para_env%ionode) THEN
         CALL open_file(file_name=scan_env%tip_cube_file, &
                        file_status="OLD", &
                        file_form="FORMATTED", &
                        file_action="READ", &
                        unit_number=extunit)
         !skip header comments
         DO i = 1, 2
            READ (extunit, *)
         ENDDO
         READ (extunit, *) nat, rdum
         DO i = 1, 3
            READ (extunit, *) npts(i), dcell(i, 1:3)
            dcell(i, 1:3) = npts(i)*dcell(i, 1:3)
         ENDDO
         CALL close_file(unit_number=extunit)
      END IF

      CALL mp_bcast(npts, para_env%source, para_env%group)
      CALL mp_bcast(dcell, para_env%source, para_env%group)

      NULLIFY (pw, pw_grid)
      CALL pw_grid_create(pw_grid, para_env%group)
      CALL pw_grid_setup(dcell, pw_grid, npts=npts)
      CALL pw_create(pw, pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
!deb
      scaling = 0.1_dp
!deb
      CALL cp_cube_to_pw(pw, scan_env%tip_cube_file, scaling, silent=.TRUE.)
      scan_env%tip_pw_r => pw
      NULLIFY (pw)
      CALL pw_create(pw, pw_grid, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_transfer(scan_env%tip_pw_r, pw)
      scan_env%tip_pw_g => pw
      CALL pw_grid_release(pw_grid)

   END SUBROUTINE read_tip_file

END MODULE tip_scan_methods
