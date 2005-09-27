!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/wiener_process [1.0] *
!!
!!   NAME
!!     wiener_process
!!
!!   FUNCTION
!!     Handling of the Wiener process currently employed in turn of the
!!     Langevin dynamics.
!!
!!   AUTHOR
!!     Matthias Krack (05.07.2005)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE wiener_process

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type
  USE input_constants,                 ONLY: langevin_ensemble,&
                                             md_all_restart,&
                                             md_init_start,&
                                             md_pos_restart,&
                                             md_pos_vel_restart
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: mdio_parameters_type,&
                                             simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             dump_rng_stream,&
                                             next_rng_seed,&
                                             read_rng_stream
  USE particle_list_types,             ONLY: particle_list_type
  USE string_utilities,                ONLY: ascii_to_string,&
                                             string_to_ascii
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

  ! Global parameters in this module

  CHARACTER(LEN=*), PARAMETER :: moduleN = "wiener_process"

  ! Public subroutines

  PUBLIC :: create_wiener_process,&
            read_md_restart_file,&
            write_md_restart_file

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE create_wiener_process(md_env,error)

    ! Purpose: Create a Wiener process for Langevin dynamics and initialize an
    !          independent random number generator for each atom in all force
    !          environment and all the subsystems/fragments therein.

    ! History: - Creation (06.07.2005,MK)

    ! *************************************************************************

    TYPE(md_environment_type), POINTER           :: md_env
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    ! Local parameters

    CHARACTER(LEN=*), PARAMETER :: routineN = "create_wiener_process",&
                                   routineP = moduleN//":"//routineN

    ! Local variables

    TYPE(atomic_kind_list_type), POINTER      :: atomic_kinds
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(distribution_1d_type), POINTER       :: local_particles
    TYPE(particle_list_type), POINTER         :: particles
    TYPE(simulation_parameters_type), POINTER :: simpar

    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(force_env_p_type), DIMENSION(:), POINTER    :: force_env_p

    REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: seed

    REAL(KIND=dp), DIMENSION(3,2) :: initial_seed

    CHARACTER(LEN=404) :: string
    CHARACTER(LEN=40)  :: name

    INTEGER :: i,iparticle,iparticle_kind,iparticle_local,istat,isubsys,&
               nparticle,nparticle_kind,nparticle_local,unit_number
    LOGICAL :: failure

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ascii

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition (ASSOCIATED(md_env),cp_failure_level,routineP,error,failure)

    CALL get_md_env(md_env=md_env,&
                    force_env_p=force_env_p,&
                    para_env=para_env,&
                    simpar=simpar)

    IF (simpar%ensemble == langevin_ensemble) THEN

      ! Load initial seed (not needed for a restart)

      initial_seed = next_rng_seed()

      DO i=1,SIZE(force_env_p)

        CALL force_env_get(force_env_p(i)%force_env,subsys=subsys)

        DO isubsys=1,SIZE(subsys)

          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             local_particles=local_particles,&
                             particles=particles,&
                             error=error)

          nparticle_kind = atomic_kinds%n_els
          nparticle = particles%n_els

          ! Allocate the (local) data structures for the Wiener process

          ALLOCATE(local_particles%local_particle_set(nparticle_kind),&
                   STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          DO iparticle_kind=1,nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            ALLOCATE (local_particles%local_particle_set(iparticle_kind)%&
                      &rng(nparticle_local),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            DO iparticle_local=1,nparticle_local
              NULLIFY (local_particles%local_particle_set(iparticle_kind)%&
                       &rng(iparticle_local)%stream)
            END DO
          END DO

          ! Each process generates all seeds. The seed generation should be
          ! quite fast and in this way a broadcast is avoided.

          ALLOCATE (seed(3,2,nparticle),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          seed(:,:,1) = initial_seed
          DO iparticle=2,nparticle
            seed(:,:,iparticle) = next_rng_seed(seed(:,:,iparticle-1))
          END DO

          ! Update initial seed

          initial_seed = next_rng_seed(seed(:,:,nparticle))

          ! Create a random number stream (Wiener process) for each particle

          DO iparticle_kind=1,nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%&
                          &array(iparticle_local)
              WRITE (UNIT=name,FMT="(A,I6)") "RNG for particle",iparticle
              CALL create_rng_stream(local_particles%&
                                     &local_particle_set(iparticle_kind)%&
                                     &rng(iparticle_local)%stream,&
                                     name=name,&
                                     extended_precision=.TRUE.,&
                                     seed=seed(:,:,iparticle))
            END DO
          END DO

          DEALLOCATE (seed,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        END DO

      END DO

    END IF

  END SUBROUTINE create_wiener_process

  ! ***************************************************************************

  SUBROUTINE read_md_restart_file(md_env,error)

    ! Purpose: Read the data needed/requested for a restart of a MD run from
    !          the restart unit unit_number.

    ! History: - Creation (11.07.2005,MK)

    ! *************************************************************************

    TYPE(md_environment_type), POINTER           :: md_env
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    ! Local parameters

    CHARACTER(LEN=*), PARAMETER :: routineN = "read_md_restart_file",&
                                   routineP = moduleN//":"//routineN

    ! Local variables

    TYPE(atomic_kind_list_type), POINTER      :: atomic_kinds
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(distribution_1d_type), POINTER       :: local_particles
    TYPE(mdio_parameters_type), POINTER       :: mdio
    TYPE(particle_list_type), POINTER         :: particles
    TYPE(simulation_parameters_type), POINTER :: simpar

    CHARACTER(LEN=404) :: string
    CHARACTER(LEN=240) :: end_section,section

    INTEGER :: i,inumber,iparticle,iparticle_kind,iparticle_local,istat,&
               isubsys,nparticle,nparticle_kind,nparticle_local,unit_number
    LOGICAL :: failure,section_read

    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(force_env_p_type), DIMENSION(:), POINTER    :: force_env_p

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ascii

    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: buffer

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition (ASSOCIATED(md_env),cp_failure_level,routineP,error,failure)

    CALL get_md_env(md_env=md_env,&
                    force_env_p=force_env_p,&
                    mdio=mdio,&
                    para_env=para_env,&
                    simpar=simpar)

    ! Quick return, if no MD restart is requested

    IF (simpar%read_type == md_init_start) RETURN

    IF (para_env%mepos == para_env%source) THEN
      CALL open_file(file_name=mdio%dump_file_name,&
                     file_action="READ",&
                     file_form="FORMATTED",&
                     file_status="OLD",&
                     unit_number=unit_number)
    END IF

    IF (ASSOCIATED(force_env_p)) THEN

      DO i=1,SIZE(force_env_p)

        IF (para_env%mepos == para_env%source) THEN
          READ (UNIT=unit_number,FMT=*) section,inumber
          IF ((TRIM(section) /= "&FORCE_ENVIRONMENT").OR.&
              (inumber /= i)) THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                           routineP,"Found """//TRIM(section)//&
                           """, but &FORCE_ENVIRONMENT was expected",error,&
                           failure)
          END IF
        END IF

        CALL force_env_get(force_env_p(i)%force_env,subsys=subsys)

        DO isubsys=1,SIZE(subsys)

          IF (para_env%mepos == para_env%source) THEN
            READ (UNIT=unit_number,FMT=*) section,inumber
            IF ((TRIM(section) /= "&SUBSYSTEM").OR.(inumber /= isubsys)) THEN
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                             routineP,"Found """//TRIM(section)//&
                             """, but &SUBSYSTEM was expected",error,&
                             failure)
            END IF
          END IF

          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             local_particles=local_particles,&
                             particles=particles,&
                             error=error)

          nparticle_kind = atomic_kinds%n_els
          nparticle = particles%n_els

          DO
            section_read = .FALSE.
            IF (para_env%mepos == para_env%source) THEN
              READ (UNIT=unit_number,FMT=*) section
              SELECT CASE (TRIM(section))
              CASE ("&COORDINATES","&VELOCITIES")
                ALLOCATE (buffer(3,nparticle),STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                DO iparticle=1,nparticle
                  READ (UNIT=unit_number,FMT=*) buffer(1:3,iparticle)
                END DO
                section_read = .TRUE.
              CASE ("&WIENER_PROCESS")
                ALLOCATE (ascii(LEN(string),nparticle),STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                ascii = 0
                DO iparticle=1,nparticle
                  READ (UNIT=unit_number,FMT="(A)") string
                  CALL string_to_ascii(string,ascii(:,iparticle))
                END DO
                section_read = .TRUE.
              CASE DEFAULT
                BACKSPACE (UNIT=unit_number)
              END SELECT
              IF (section_read) THEN
                READ (UNIT=unit_number,FMT=*) end_section
                IF (TRIM(end_section) /= "&END") THEN
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                                 routineP,"Found """//TRIM(end_section)//&
                                 """, but &END was expected",error,&
                                 failure)
                END IF
              END IF
            END IF

            CALL mp_bcast(section,para_env%source,para_env%group)

            SELECT CASE (TRIM(section))
            CASE ("&COORDINATES")
              SELECT CASE (simpar%read_type)
              CASE (md_pos_restart,md_pos_vel_restart,md_all_restart)
                IF (para_env%mepos /= para_env%source) THEN
                  ALLOCATE (buffer(3,nparticle),STAT=istat)
                  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                END IF
                CALL mp_bcast(buffer,para_env%source,para_env%group)
                DO iparticle=1,nparticle
                  particles%els(iparticle)%r(1:3) = buffer(1:3,iparticle)
                END DO
                DEALLOCATE (buffer,STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              END SELECT
            CASE ("&VELOCITIES")
              SELECT CASE (simpar%read_type)
              CASE (md_pos_vel_restart,md_all_restart)
                IF (para_env%mepos /= para_env%source) THEN
                  ALLOCATE (buffer(3,nparticle),STAT=istat)
                  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                END IF
                CALL mp_bcast(buffer,para_env%source,para_env%group)
                DO iparticle=1,nparticle
                  particles%els(iparticle)%v(1:3) = buffer(1:3,iparticle)
                END DO
                DEALLOCATE (buffer,STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              END SELECT
            CASE ("&WIENER_PROCESS")
              SELECT CASE (simpar%read_type)
              CASE (md_all_restart)
                IF (para_env%mepos /= para_env%source) THEN
                  ALLOCATE (ascii(LEN(string),nparticle),STAT=istat)
                  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                  ascii = 0
                END IF
                CALL mp_bcast(ascii,para_env%source,para_env%group)
                DO iparticle_kind=1,nparticle_kind
                  nparticle_local = local_particles%n_el(iparticle_kind)
                  DO iparticle_local=1,nparticle_local
                    iparticle = local_particles%list(iparticle_kind)%&
                                &array(iparticle_local)
                    IF (iparticle == local_particles%list(iparticle_kind)%&
                                     &array(iparticle_local)) THEN
                      CALL ascii_to_string(ascii(:,iparticle),string)
                      CALL read_rng_stream(rng_stream=local_particles%&
                                           &local_particle_set(iparticle_kind)%&
                                           &rng(iparticle_local)%stream,&
                                           rng_record=string,&
                                           error=error)
                    END IF
                  END DO
                END DO
                DEALLOCATE (ascii,STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              END SELECT
            CASE DEFAULT
              EXIT
            END SELECT

          END DO

          IF (para_env%mepos == para_env%source) THEN
            READ (UNIT=unit_number,FMT=*) end_section
            IF (TRIM(end_section) /= "&END") THEN
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                             routineP,"Found """//TRIM(end_section)//&
                             """, but &END was expected",error,&
                             failure)
            END IF
          END IF

        END DO

        IF (para_env%mepos == para_env%source) THEN
          READ (UNIT=unit_number,FMT=*) end_section
          IF (TRIM(end_section) /= "&END") THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                           routineP,"Found """//TRIM(end_section)//&
                           """, but &END was expected",error,&
                           failure)
          END IF
        END IF

      END DO

    END IF

    IF (para_env%mepos == para_env%source) THEN
      CALL close_file(unit_number=unit_number)
    END IF

  END SUBROUTINE read_md_restart_file

  ! ***************************************************************************

  SUBROUTINE write_md_restart_file(md_env,hmc,error)

    ! Purpose: Write the data needed for a restart of a MD run to the restart
    !          unit unit_number.

    ! History: - Creation (01.07.2005,MK)

    ! *************************************************************************

    TYPE(md_environment_type), POINTER           :: md_env
    LOGICAL, OPTIONAL                            :: hmc
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    ! Local parameters

    CHARACTER(LEN=*), PARAMETER :: routineN = "write_md_restart_file",&
                                   routineP = moduleN//":"//routineN

    ! Local variables

    TYPE(atomic_kind_list_type), POINTER      :: atomic_kinds
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(distribution_1d_type), POINTER       :: local_particles
    TYPE(mdio_parameters_type), POINTER       :: mdio
    TYPE(particle_list_type), POINTER         :: particles
    TYPE(simulation_parameters_type), POINTER :: simpar

    CHARACTER(LEN=404) :: string

    INTEGER :: i,iparticle,iparticle_kind,iparticle_local,istat,isubsys,&
               nparticle,nparticle_kind,nparticle_local,unit_number
    LOGICAL :: failure, hmc_flag

    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(force_env_p_type), DIMENSION(:), POINTER    :: force_env_p

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ascii

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition (ASSOCIATED(md_env),cp_failure_level,routineP,error,failure)

    IF (PRESENT(hmc)) THEN
      hmc_flag = hmc
    ELSE
      hmc_flag = .FALSE.
    END IF

    CALL get_md_env(md_env=md_env,&
                    force_env_p=force_env_p,&
                    mdio=mdio,&
                    para_env=para_env,&
                    simpar=simpar)

    IF (para_env%mepos == para_env%source) THEN
      CALL open_file(file_name=mdio%dump_file_name,&
                     file_action="WRITE",&
                     file_form="FORMATTED",&
                     file_status="REPLACE",&
                     unit_number=unit_number)
    END IF

    IF (ASSOCIATED(force_env_p)) THEN

      DO i=1,SIZE(force_env_p)

        IF (para_env%mepos == para_env%source) THEN
          WRITE (UNIT=unit_number,FMT="(A,I5)") "&FORCE_ENVIRONMENT",i
        END IF

        CALL force_env_get(force_env_p(i)%force_env,subsys=subsys)

        DO isubsys=1,SIZE(subsys)

          IF (para_env%mepos == para_env%source) THEN
            WRITE (UNIT=unit_number,FMT="(A,I5)") " &SUBSYSTEM",isubsys
          END IF

          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             local_particles=local_particles,&
                             particles=particles,&
                             error=error)

          nparticle_kind = atomic_kinds%n_els
          nparticle = particles%n_els

          IF (para_env%mepos == para_env%source) THEN

            ! Write the current particle coordinates to the restart unit

            WRITE (UNIT=unit_number,FMT="(A)") "  &COORDINATES"
            WRITE (UNIT=unit_number,FMT="(3ES25.14)")&
              (particles%els(iparticle)%r(1:3),iparticle=1,nparticle)
            WRITE (UNIT=unit_number,FMT="(A)") "  &END COORDINATES"

            ! Write the current particle velocities to the restart unit

            WRITE (UNIT=unit_number,FMT="(A)") "  &VELOCITIES"
            WRITE (UNIT=unit_number,FMT="(3ES25.14)")&
              (particles%els(iparticle)%v(1:3),iparticle=1,nparticle)
            WRITE (UNIT=unit_number,FMT="(A)") "  &END VELOCITIES"

          END IF

          IF (.NOT. hmc_flag) THEN

          ! Write current status of Wiener process to the restart unit

          ALLOCATE (ascii(LEN(string),nparticle),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ascii = 0

          DO iparticle=1,nparticle
            DO iparticle_kind=1,nparticle_kind
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                IF (iparticle == local_particles%list(iparticle_kind)%&
                                 &array(iparticle_local)) THEN
                  CALL dump_rng_stream(local_particles%&
                                       &local_particle_set(iparticle_kind)%&
                                       rng(iparticle_local)%stream,string)
                  CALL string_to_ascii(string,ascii(:,iparticle))
                END IF
              END DO
            END DO
          END DO

          CALL mp_sum(ascii,para_env%group)

          IF (para_env%mepos == para_env%source) THEN
            WRITE (UNIT=unit_number,FMT="(A)") "  &WIENER_PROCESS"
            DO iparticle=1,nparticle
              CALL ascii_to_string(ascii(:,iparticle),string)
              WRITE (UNIT=unit_number,FMT="(A)") TRIM(string)
            END DO
            WRITE (UNIT=unit_number,FMT="(A)") "  &END WIENER_PROCESS"
          END IF

          DEALLOCATE (ascii,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          END IF

          IF (para_env%mepos == para_env%source) THEN
            WRITE (UNIT=unit_number,FMT="(A,I5)") " &END SUBSYSTEM",isubsys
          END IF

        END DO

        IF (para_env%mepos == para_env%source) THEN
          WRITE (UNIT=unit_number,FMT="(A,I5)") "&END FORCE_ENVIRONMENT",i
        END IF

      END DO

    END IF

    IF (para_env%mepos == para_env%source) THEN
      CALL close_file(unit_number=unit_number)
    END IF

  END SUBROUTINE write_md_restart_file

 ! ***************************************************************************

END MODULE wiener_process
