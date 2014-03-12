!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
MODULE csvr_system_mapping

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE csvr_system_types,               ONLY: csvr_system_type,&
                                             csvr_thermo_create
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: debug_isotropic_limit,&
                                             map_info_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_thermo_communication, do_thermo_no_communication, &
       do_thermo_only_master, isokin_ensemble, langevin_ensemble, &
       npe_f_ensemble, npe_i_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, nve_ensemble, &
       nvt_ensemble, reftraj_ensemble
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE thermostat_mapping,              ONLY: init_baro_map_info,&
                                             thermostat_mapping_region
  USE thermostat_types,                ONLY: thermostat_info_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_mapping'

  PUBLIC :: csvr_to_particle_mapping, csvr_to_barostat_mapping,&
            csvr_to_shell_mapping

CONTAINS

! *****************************************************************************
!> \brief Creates the thermostatting for the barostat
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE csvr_to_barostat_mapping ( simpar, csvr, error )
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_barostat_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ndeg
    LOGICAL                                  :: failure
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    SELECT CASE (simpar%ensemble)
    CASE  DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,'Never reach this point!')
    CASE ( npt_i_ensemble, npt_f_ensemble )
       map_info => csvr%map_info
       map_info%dis_type = do_thermo_only_master

       ! Counting the total number of thermostats ( 1 for both NPT_I and NPT_F )
       csvr%loc_num_csvr  = 1
       csvr%glob_num_csvr = 1
       IF ( simpar % ensemble == npt_f_ensemble ) THEN
          ndeg = 9
       ELSE
          ndeg = 1
       ENDIF

       CALL init_baro_map_info(map_info, ndeg, csvr%loc_num_csvr, error)
       CALL csvr_thermo_create(csvr, error)

       ! Now that we know how many there are stick this into csvr%nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, csvr%loc_num_csvr
          csvr%nvt(i)%nkt =  simpar%temp_baro_ext * ndeg
          csvr%nvt(i)%degrees_of_freedom = ndeg
          IF (debug_isotropic_limit) THEN
             csvr%nvt(i)%nkt =  simpar%temp_baro_ext
             csvr%nvt(i)%degrees_of_freedom = 1
          END IF
       END DO
    END SELECT

  END SUBROUTINE csvr_to_barostat_mapping

! *****************************************************************************
!> \brief Creates the thermostatting maps
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE csvr_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, csvr, para_env, gci, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_particle_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap, j, natoms_local, &
                                                stat, sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fac
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    NULLIFY (massive_atom_list, deg_of_freedom)
    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,'Unknown ensemble!')
    CASE ( nve_ensemble, isokin_ensemble, npe_f_ensemble, npe_i_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program(routineN,moduleN,__LINE__,'Never reach this point!')
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble)

       CALL setup_csvr_thermostat(csvr, thermostat_info, deg_of_freedom,&
            massive_atom_list, molecule_kind_set, local_molecules, molecule_set,&
            para_env, natoms_local, simpar, sum_of_thermostats, gci, error=error)

       ! Sum up the number of degrees of freedom on each thermostat.
       ! first: initialize the target
       map_info => csvr%map_info
       map_info%s_kin = 0.0_dp
       DO i = 1, 3
          DO j = 1, natoms_local
             map_info%p_kin(i,j)%point = map_info%p_kin(i,j)%point + 1
          END DO
       END DO

       ! If thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF (map_info%dis_type == do_thermo_communication ) CALL mp_sum (map_info%s_kin, para_env%group )

       ! We know the total number of system thermostats.
       IF ( (sum_of_thermostats == 1).AND.(map_info%dis_type /= do_thermo_no_communication)) THEN
          fac = map_info%s_kin(1)-deg_of_freedom(1)-simpar%nfree_rot_transl
          IF ( fac == 0.0_dp ) THEN
             CALL stop_program(routineN,moduleN,__LINE__,'Zero degrees of freedom. Nothing to thermalize!')
          END IF
          csvr%nvt(1)%nkt = simpar%temp_ext*fac
          csvr%nvt(1)%degrees_of_freedom = FLOOR(fac)
       ELSE
          DO i = 1, csvr%loc_num_csvr
             imap = map_info%map_index(i)
             fac = (map_info%s_kin(imap) - deg_of_freedom(i) )
             csvr%nvt(i)%nkt = simpar%temp_ext * fac
             csvr%nvt(i)%degrees_of_freedom = FLOOR(fac)
          END DO
       END IF

       DEALLOCATE (deg_of_freedom,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( massive_atom_list, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE csvr_to_particle_mapping

! *****************************************************************************
!> \brief Main general setup for CSVR thermostats
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
  SUBROUTINE setup_csvr_thermostat ( csvr, thermostat_info, deg_of_freedom,&
       massive_atom_list, molecule_kind_set, local_molecules, molecule_set,&
       para_env, natoms_local, simpar, sum_of_thermostats, gci, shell, error)

    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(OUT)                     :: natoms_local
    TYPE(simpar_type), POINTER               :: simpar
    INTEGER, INTENT(OUT)                     :: sum_of_thermostats
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_csvr_thermostat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nkind, number, region
    LOGICAL                                  :: do_shell, failure
    TYPE(map_info_type), POINTER             :: map_info

    failure  = .FALSE.
    do_shell = .FALSE.
    IF (PRESENT(shell)) do_shell = shell
    map_info    => csvr%map_info

    nkind = SIZE(molecule_kind_set)
    sum_of_thermostats = thermostat_info%sum_of_thermostats
    map_info%dis_type  = thermostat_info%dis_type
    number             = thermostat_info%number_of_thermostats
    region             = csvr%region

    CALL thermostat_mapping_region(map_info, deg_of_freedom, massive_atom_list,&
       molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local,&
       simpar, number, region, gci, do_shell, thermostat_info%map_loc_thermo_gen,&
       sum_of_thermostats, error)

    ! This is the local number of available thermostats
    csvr%loc_num_csvr  = number
    csvr%glob_num_csvr = sum_of_thermostats
    CALL csvr_thermo_create(csvr, error=error)

  END SUBROUTINE setup_csvr_thermostat

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
  SUBROUTINE csvr_to_shell_mapping( thermostat_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, csvr, para_env, gci, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_shell_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap, j, nshell_local, &
                                                stat, sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_shell_list
    LOGICAL                                  :: failure
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    NULLIFY (massive_shell_list, deg_of_freedom)

    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,'Unknown ensemble!')
    CASE ( isokin_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program(routineN,moduleN,__LINE__,'Never reach this point!')
    CASE ( nve_ensemble, npe_f_ensemble, npe_i_ensemble, nvt_ensemble, npt_i_ensemble, npt_f_ensemble )

       CALL setup_csvr_thermostat(csvr, thermostat_info, deg_of_freedom, massive_shell_list,&
            molecule_kind_set, local_molecules, molecule_set, para_env, nshell_local, &
            simpar, sum_of_thermostats, gci, shell=.TRUE., error=error)

       map_info => csvr%map_info
       ! Sum up the number of degrees of freedom on each thermostat.
       ! first: initialize the target
       map_info%s_kin = 0.0_dp
       DO j = 1, nshell_local
          DO i = 1, 3
             map_info%p_kin ( i, j ) % point = map_info%p_kin ( i, j ) % point + 1
          END DO
       END DO

       ! If thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF ( map_info%dis_type == do_thermo_communication ) CALL mp_sum (map_info%s_kin, para_env%group )

       ! Now that we know how many there are stick this into csvr%nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, csvr%loc_num_csvr
          imap = map_info%map_index(i)
          csvr%nvt(i)%nkt = simpar%temp_sh_ext * map_info%s_kin(imap)
          csvr%nvt(i)%degrees_of_freedom = FLOOR(map_info%s_kin(imap))
       END DO

       DEALLOCATE (deg_of_freedom,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( massive_shell_list, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE csvr_to_shell_mapping

END MODULE csvr_system_mapping
