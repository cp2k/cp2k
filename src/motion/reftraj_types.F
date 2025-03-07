!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief initialization of the reftraj structure used to analyse
!>     previously generated trajectories
!> \par History
!>      Created 10-07 [MI]
!> \author MI
! **************************************************************************************************
MODULE reftraj_types

   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: reftraj_type, reftraj_msd_type, &
             create_reftraj, release_reftraj

   INTEGER, PARAMETER, PUBLIC :: REFTRAJ_EVAL_NONE = 101
   INTEGER, PARAMETER, PUBLIC :: REFTRAJ_EVAL_ENERGY = 102
   INTEGER, PARAMETER, PUBLIC :: REFTRAJ_EVAL_ENERGY_FORCES = 103

! **************************************************************************************************
!> \brief parameters related to the analysis of previously generated trajecorties
!> \author MI
! **************************************************************************************************
   TYPE reftraj_info_type
      INTEGER                                  :: first_snapshot = 0
      INTEGER                                  :: last_snapshot = 0
      INTEGER                                  :: stride = 0
      INTEGER                                  :: eval = REFTRAJ_EVAL_NONE
      LOGICAL                                  :: variable_volume = .FALSE.
      LOGICAL                                  :: msd = .FALSE.
      TYPE(cp_parser_type), POINTER            :: traj_parser => NULL()
      TYPE(cp_parser_type), POINTER            :: cell_parser => NULL()
   END TYPE reftraj_info_type

! **************************************************************************************************
   TYPE reftraj_msd_type
   LOGICAL                                  :: disp_atom = .FALSE., msd_kind = .FALSE., msd_molecule = .FALSE., msd_region = .FALSE.
      INTEGER                                  :: num_disp_atom = 0, ref0_unit = 0
      INTEGER, POINTER, DIMENSION(:)           :: disp_atom_index => NULL()
    REAL(KIND=dp)                            :: disp_atom_tol = 0.0_dp, drcom(3) = 0.0_dp, ref0_com(3) = 0.0_dp, total_mass = 0.0_dp
      REAL(KIND=dp), POINTER, DIMENSION(:, :)   :: disp_atom_dr => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)   :: ref0_pos => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)   :: ref0_com_molecule => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)   :: val_msd_kind => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)   :: val_msd_molecule => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)   :: val_msd_region => NULL()
   END TYPE reftraj_msd_type

! **************************************************************************************************
   TYPE reftraj_type
      INTEGER                                  :: itimes = 0
      INTEGER                                  :: itimes0 = 0
      INTEGER                                  :: isnap = 0
      INTEGER                                  :: natom = 0
      LOGICAL                                  :: init = .FALSE.
      REAL(KIND=dp)                            :: epot = 0.0_dp, epot0 = 0.0_dp, time = 0.0_dp, time0 = 0.0_dp
      TYPE(reftraj_info_type), POINTER         :: info => NULL()
      TYPE(reftraj_msd_type), POINTER          :: msd => NULL()
   END TYPE reftraj_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reftraj_types'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param reftraj ...
!> \param reftraj_section ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE create_reftraj(reftraj, reftraj_section, para_env)

      TYPE(reftraj_type), INTENT(OUT)                    :: reftraj
      TYPE(section_vals_type), POINTER                   :: reftraj_section
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=default_path_length)                 :: filename
      LOGICAL                                            :: old_eval_ef, old_eval_forces

      NULLIFY (reftraj%info)
      NULLIFY (reftraj%msd)

      ALLOCATE (reftraj%info)
      NULLIFY (reftraj%info%traj_parser)
      NULLIFY (reftraj%info%cell_parser)

      ! Initialize parser for trajectory
      CALL section_vals_val_get(reftraj_section, "TRAJ_FILE_NAME", c_val=filename)
      ALLOCATE (reftraj%info%traj_parser)
      CALL parser_create(reftraj%info%traj_parser, filename, para_env=para_env)

      CALL section_vals_val_get(reftraj_section, "VARIABLE_VOLUME", l_val=reftraj%info%variable_volume)
      IF (reftraj%info%variable_volume) THEN
         ! In case requested initialize parser for cell
         CALL section_vals_val_get(reftraj_section, "CELL_FILE_NAME", c_val=filename)
         ALLOCATE (reftraj%info%cell_parser)
         CALL parser_create(reftraj%info%cell_parser, filename, para_env=para_env)
      END IF

      CALL section_vals_val_get(reftraj_section, "FIRST_SNAPSHOT", i_val=reftraj%info%first_snapshot)
      CALL section_vals_val_get(reftraj_section, "LAST_SNAPSHOT", i_val=reftraj%info%last_snapshot)
      CALL section_vals_val_get(reftraj_section, "STRIDE", i_val=reftraj%info%stride)
      CALL section_vals_val_get(reftraj_section, "EVAL", i_val=reftraj%info%eval)

      ! Read deprecated keywords to retain backwards compatibility.
      ! For details see: https://github.com/cp2k/cp2k/issues/894
      CALL section_vals_val_get(reftraj_section, "EVAL_ENERGY_FORCES", l_val=old_eval_ef)
      CALL section_vals_val_get(reftraj_section, "EVAL_FORCES", l_val=old_eval_forces)
      IF (old_eval_ef) reftraj%info%eval = REFTRAJ_EVAL_ENERGY
      IF (old_eval_forces) reftraj%info%eval = REFTRAJ_EVAL_ENERGY_FORCES

      CALL section_vals_val_get(reftraj_section, "MSD%_SECTION_PARAMETERS_", &
                                l_val=reftraj%info%msd)

   END SUBROUTINE create_reftraj

! **************************************************************************************************
!> \brief ...
!> \param reftraj ...
!> \par History
!>      10.2007 created
!> \author MI
! **************************************************************************************************
   SUBROUTINE release_reftraj(reftraj)

      TYPE(reftraj_type), INTENT(INOUT)                  :: reftraj

      IF (ASSOCIATED(reftraj%info%traj_parser)) THEN
         CALL parser_release(reftraj%info%traj_parser)
         DEALLOCATE (reftraj%info%traj_parser)
      END IF
      IF (ASSOCIATED(reftraj%info%cell_parser)) THEN
         CALL parser_release(reftraj%info%cell_parser)
         DEALLOCATE (reftraj%info%cell_parser)
      END IF
      IF (ASSOCIATED(reftraj%info)) THEN
         DEALLOCATE (reftraj%info)
      END IF
      IF (ASSOCIATED(reftraj%msd)) THEN
         DEALLOCATE (reftraj%msd%ref0_pos)
         IF (reftraj%msd%msd_kind) THEN
            DEALLOCATE (reftraj%msd%val_msd_kind)
         END IF
         IF (reftraj%msd%msd_molecule) THEN
            DEALLOCATE (reftraj%msd%val_msd_molecule)
            DEALLOCATE (reftraj%msd%ref0_com_molecule)
         END IF
         IF (reftraj%msd%disp_atom) THEN
            DEALLOCATE (reftraj%msd%disp_atom_index)
            DEALLOCATE (reftraj%msd%disp_atom_dr)
         END IF

         DEALLOCATE (reftraj%msd)
      END IF

   END SUBROUTINE release_reftraj

END MODULE reftraj_types
