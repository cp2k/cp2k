!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/csvr_system_mapping [1.0] *
!!
!!   NAME
!!     csvr_system_mapping
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] 10.2007- University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************

MODULE csvr_system_mapping

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE csvr_system_types,               ONLY: csvr_system_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             get_dyn_coeff
  USE extended_system_types,           ONLY: debug_isotropic_limit,&
                                             map_info_type
  USE input_constants,                 ONLY: &
       do_region_global, do_region_massive, do_thermo_communication, &
       do_thermo_no_communication, do_thermo_only_master, isokin_ensemble, &
       langevin_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, nve_ensemble, &
       nvt_ensemble, reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE termination,                     ONLY: stop_program
  USE thermostat_mapping,              ONLY: init_baro_map_info,&
                                             init_coeff_map_info,&
                                             thermostat_mapping_region
  USE thermostat_types,                ONLY: thermostat_info_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_mapping'

  PUBLIC :: csvr_to_particle_mapping, csvr_to_barostat_mapping,&
            csvr_to_coef_mapping, csvr_to_shell_mapping

CONTAINS
!******************************************************************************
!!****** csvr_system_mapping/csvr_to_barostat_mapping [1.0] *
!!
!!   NAME
!!     csvr_to_barostat_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting for the barostat
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] 10.2007- University of Zurich
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** *************************************************************************
  SUBROUTINE csvr_to_barostat_mapping ( simpar, csvr, para_env, error )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_barostat_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, number, stat
    LOGICAL                                  :: failure
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    SELECT CASE (simpar%ensemble)
    CASE  DEFAULT
       CALL stop_program(routineN, 'Never reach this point!' )
    CASE ( npt_i_ensemble, npt_f_ensemble )
       map_info => csvr%map_info
       map_info%dis_type = do_thermo_only_master

       ! Counting the total number of thermostats ( 1 for both NPT_I and NPT_F )
       csvr%num_csvr = 1
       IF ( simpar % ensemble == npt_f_ensemble ) THEN
          number = 9
       ELSE
          number = 1
       ENDIF

       CALL init_baro_map_info(map_info, number, csvr%num_csvr, error)
 
       ALLOCATE ( csvr%nvt(csvr%num_csvr),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Now that we know how many there are stick this into csvr%nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, csvr%num_csvr
          csvr%nvt(i)%nkt = simpar % temp_ext * number
          IF (debug_isotropic_limit) THEN
             csvr%nvt(i)%nkt = simpar % temp_ext
          END IF
       END DO
    END SELECT
    
  END SUBROUTINE csvr_to_barostat_mapping

!******************************************************************************
!!****** csvr_system_mapping/csvr_to_particle_mapping [1.0] *
!!
!!   NAME
!!     csvr_to_particle_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting maps
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] 10.2007- University of Zurich
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** *************************************************************************
  SUBROUTINE csvr_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, csvr, para_env, gci, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_particle_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, natoms_local, stat, &
                                                sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fac
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    NULLIFY (massive_atom_list, deg_of_freedom)
    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program ( routineN, 'Unknown ensemble!' )
    CASE ( nve_ensemble, isokin_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program ( routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble )

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
             CALL stop_program ( routineN,'Zero degrees of freedom. Nothing to thermalize!' )
          END IF
          csvr%nvt(1)%nkt = simpar%temp_ext*fac
       ELSE
          DO i = 1, csvr%num_csvr
             fac = (map_info%s_kin(i) - deg_of_freedom(i) )
             csvr%nvt(i)%nkt = simpar%temp_ext * fac 
          END DO
       END IF
       
       DEALLOCATE (deg_of_freedom,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( massive_atom_list, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE csvr_to_particle_mapping

!******************************************************************************
!!****** csvr_system_mapping/setup_csvr_thermostat [1.0] *
!!
!!   NAME
!!     setup_csvr_thermostat
!!
!!   FUNCTION
!!     Main general setup for Nose-Hoover thermostats
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
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
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    INTEGER, INTENT(OUT)                     :: sum_of_thermostats
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_csvr_thermostat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nkind, number, region, stat
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
       simpar, number, region, gci, do_shell, error)

    ALLOCATE (csvr%nvt(number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Now that we know how many there are stick this into csvr%nkt
    ! (number of degrees of freedom times k_B T for the first thermostat
    !  on the chain)
    csvr % num_csvr = number

  END SUBROUTINE setup_csvr_thermostat

!******************************************************************************
!!****** csvr_system_mapping/csvr_to_coef_mapping [1.0] *
!!
!!   NAME
!!     csvr_to_coef_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting for the
!!     polarization  coefficients'  degrees of freedom
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE csvr_to_coef_mapping ( dyn_coeff_set, csvr, kg_coeff_section, para_env, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(section_vals_type), POINTER         :: kg_coeff_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_coef_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ikind, j, kg_coeff_ens, &
                                                ncoef_local, ncoef_of_kind, &
                                                nkind, number, stat
    LOGICAL                                  :: failure, nointer
    REAL(KIND=dp)                            :: kg_coeff_temp
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(section_vals_type), POINTER         :: nose_section

    failure = .FALSE.
    nose_section => section_vals_get_subs_vals(kg_coeff_section,"thermostat%nose",error=error)

    CALL section_vals_val_get(kg_coeff_section,"ensemble",i_val=kg_coeff_ens,error=error)
    CALL section_vals_val_get(kg_coeff_section,"temperature",r_val=kg_coeff_temp,error=error)

    SELECT CASE ( kg_coeff_ens )
    CASE  DEFAULT
       CALL stop_program(routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble )
       number      = 0
       ncoef_local = 0
       nointer     = .TRUE.
       nkind       = SIZE ( dyn_coeff_set%coeffs_of_kind)
       map_info    => csvr%map_info
       DO ikind=1,nkind
          IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
             CALL get_dyn_coeff(coeffs=dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                  ncoef_of_kind=ncoef_of_kind,error=error)
             ncoef_local= ncoef_local + ncoef_of_kind
          ENDIF
       END DO
       
       ! Counting the global number of thermostats
       IF ( csvr%region == do_region_global )  THEN
          nointer=.FALSE.
          number=1
       ELSEIF ( csvr%region == do_region_massive )  THEN
          number=ncoef_local
       ELSE
          CALL stop_program ( routineN, "Unknown region for KG COEFF!" )
       ENDIF

       ! Coefficients are local. communications needed only with GLOBAL
       IF ( nointer ) THEN
          map_info%dis_type = do_thermo_no_communication
       ELSE
          map_info%dis_type = do_thermo_communication
       END IF
       csvr % num_csvr = number
       CALL init_coeff_map_info(map_info, number, csvr%region, ncoef_local, error)

       ALLOCATE ( csvr%nvt(csvr%num_csvr),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Find the number of degrees of freedom on each thermostat
       map_info%s_kin = 0.0_dp
       DO i = 1, ncoef_local
          DO j = 1, SIZE ( map_info%p_kin, 1 )
             map_info%p_kin ( j, i ) % point  = map_info%p_kin ( j, i ) % point  + 1.0_dp
          END DO
       ENDDO

       ! if thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF ( map_info%dis_type == do_thermo_communication ) CALL mp_sum (map_info%s_kin, para_env%group )

       !
       ! Now that we know how many there are stick this into csvr%nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, csvr% num_csvr
          csvr%nvt(i)%nkt = kg_coeff_temp * map_info%s_kin ( i )
       END DO

    END SELECT

  END SUBROUTINE csvr_to_coef_mapping

!******************************************************************************
!!****** csvr_system_mapping/csvr_to_shell_mapping [1.0] *
!!
!!   NAME
!!     csvr_to_shell_mapping
!!
!!   FUNCTION
!!     
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!   
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE csvr_to_shell_mapping( thermostat_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, csvr, para_env, gci, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_to_shell_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nshell_local, stat, &
                                                sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_shell_list
    LOGICAL                                  :: failure
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    NULLIFY (massive_shell_list, deg_of_freedom)

    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program ( routineN, 'Unknown ensemble!' )
    CASE ( nve_ensemble, isokin_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program ( routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble )

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
       DO i = 1, csvr%num_csvr
          csvr%nvt(i)%nkt = simpar%temp_sh_ext * map_info%s_kin ( i )
       END DO

       DEALLOCATE (deg_of_freedom,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( massive_shell_list, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE csvr_to_shell_mapping

END MODULE csvr_system_mapping
