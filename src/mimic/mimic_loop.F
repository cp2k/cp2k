!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module containing the main loop for a MiMiC run 
!> \par History
!>      03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_loop

   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE force_env_types,                 ONLY: force_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mimic_communicator,              ONLY: mimic_communicator_type

#if defined(__MIMIC)
   USE mcl
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   
   PRIVATE

   PUBLIC :: do_mimic_loop

CONTAINS

! **************************************************************************************************
!> \brief The main loop for a MiMiC run
!> \param force_env ...
!> \par History
!>      3.2023 created [AA]
! **************************************************************************************************
   SUBROUTINE do_mimic_loop(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'do_mimic_loop'

      TYPE(mimic_communicator_type)                      :: mimic_comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: is_last_step, executed
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(force_env))
      
      NULLIFY (para_env)
      
      para_env => force_env%para_env

#if defined(__MIMIC)
      CALL mimic_comm%initialize(force_env)

      is_last_step = .FALSE.
      DO WHILE (.NOT. is_last_step)
         !> Receive command from MiMiC server
         CALL mimic_comm%receive_command()
         
         !> Execute the received command (communication)
         CALL mimic_comm%execute_command(executed)
            
         !> Execute the received command (computation)
         IF (.NOT. executed) THEN
            IF (mimic_comm%is_cmd(MCL_COMPUTE_FORCES)) THEN
               CALL force_env_calc_energy_force(force_env, calc_force=.TRUE.)
            ELSE IF (mimic_comm%is_cmd(MCL_EXIT)) THEN
               is_last_step = .TRUE.
            ELSE
              CPABORT("Unrecognized MiMiC request!")
            END IF
         END IF
      END DO
      CALL para_env%sync()

      CALL mimic_comm%finalize()

#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

      CALL timestop(handle)

   END SUBROUTINE do_mimic_loop

END MODULE mimic_loop

