!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/nhc_mapping [1.0] *
!!
!!   NAME
!!     nhc
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM, 20-Feb-01
!!     JGH (10-Mar-2001)
!!     CJM (10-Apr-2001)
!!
!!   SOURCE
!******************************************************************************

MODULE nhc_mapping

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             get_dyn_coeff
  USE extended_system_types,           ONLY: lnhc_parameters_type
  USE init_extended_system_variables,  ONLY: debug_isotropic_limit
  USE input_constants,                 ONLY: &
       do_region_global, do_region_massive, do_region_molecule, &
       do_thermo_communication, do_thermo_no_communication, &
       do_thermo_only_master, isokin_ensemble, langevin_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_ensemble, reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE termination,                     ONLY: stop_program
  USE thermostat_mapping,              ONLY: massive_list_generate,&
                                             thermostat_mapping_region
  USE thermostat_types,                ONLY: thermostats_info_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nhc_mapping'

  PUBLIC :: nhc_to_particle_mapping, nhc_to_barostat_mapping, &
       nhc_to_coef_mapping, nhc_to_shell_mapping

CONTAINS
!******************************************************************************
!!****** nhc_mapping/nhc_to_barostat_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_barostat_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting for the barostat
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     CJM, 20-Feb-01  : nhc structure allocated to zero when not in use
!!     JGH (10-Mar-2001) : set nhc variables to zero when not in use
!!
!!*** *************************************************************************
  SUBROUTINE nhc_to_barostat_mapping ( simpar, nhc, para_env, error )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_barostat_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, number, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! begin local definitions
    nhc % dis_type = do_thermo_only_master
    SELECT CASE (simpar%ensemble)
    CASE  DEFAULT
       CALL stop_program(routineN, 'Never reach this point!' )
    CASE ( npt_i_ensemble, npt_f_ensemble )
       nhc % nyosh = simpar % nyosh
       nhc % nc = simpar % nc
       nhc % nhc_len = simpar % nhclen
       
       ! Counting the total number of thermostats ( 1 for both NPT_I and NPT_F )
       nhc % num_nhc = 1
       IF ( simpar % ensemble == npt_f_ensemble ) THEN
          number = 9
       ELSE
          number = 1
       ENDIF

       ! Now we can allocate the target array s_kin and p_kin
       ALLOCATE ( nhc % s_kin ( nhc % num_nhc ), STAT = stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % v_scale ( nhc % num_nhc ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % p_kin ( 1, number ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % p_scale ( 1, number ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Allocate the index array
       ALLOCATE ( nhc%index(1), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Begin the mapping loop
       DO i = 1, number
          nhc % p_kin(1,i)%point   => nhc%s_kin(1)
          nhc % p_scale(1,i)%point => nhc%v_scale(1)
       END DO
       nhc%index(1) = 1

       ALLOCATE ( nhc%nvt(nhc%nhc_len, nhc%num_nhc), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Now that we know how many there are stick this into nhc % static % nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, nhc % num_nhc
          nhc % nvt ( 1, i ) % nkt = simpar % temp_ext * number
          IF (debug_isotropic_limit) THEN
             nhc % nvt ( 1, i ) % nkt = simpar % temp_ext
          END IF
       END DO

       ! getting the number of degrees of freedom times k_B T for the rest of the chain
       DO i = 2, nhc % nhc_len
          nhc % nvt ( i, : ) %nkt = simpar % temp_ext
       END DO
    END SELECT

  END SUBROUTINE nhc_to_barostat_mapping

!******************************************************************************
!!****** nhc_mapping/nhc_to_particle_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_particle_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting maps
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     29-Nov-00 (JGH) correct counting of DOF if constraints are off
!!     CJM, 20-Feb-01  : nhc structure allocated to zero when not in use
!!     JGH (10-Mar-2001) : set nhc variables to zero when not in use
!!     CJM(10-NOV-2001) : New parallelization with new molecule structures
!!     Teodoro Laino 09.2007 [tlaino] - University of Zurich - cleaning and updating
!!
!!*** *************************************************************************
  SUBROUTINE nhc_to_particle_mapping ( thermostats_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, nhc, para_env, error)

    TYPE(thermostats_info_type), POINTER     :: thermostats_info
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_particle_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_increment, first_atom, i, icount, ielement, ii, iii, &
      ikind, imol, imol_local, ipart, j, jj, kk, last_atom, natom, &
      natoms_local, nc, nkind, nmol_local, number, region, stat, &
      sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: const_mol, deg_of_freedom, &
                                                massive_atom_list, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    LOGICAL                                  :: failure, have_global
    REAL(KIND=dp)                            :: fac
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    NULLIFY (massive_atom_list, molecule_kind, const_mol, deg_of_freedom, tot_const,&
             point)
    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program ( routineN, 'Unnknown ensemble!' )
    CASE ( nve_ensemble, isokin_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program ( routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble )
       nhc%nyosh = simpar%nyosh
       nhc%nc = simpar%nc
       nhc%nhc_len = simpar%nhclen

       nkind = SIZE ( molecule_kind_set )
       sum_of_thermostats = thermostats_info%sum_of_thermostats
       have_global        = thermostats_info%have_global
       nhc%dis_type       = thermostats_info%dis_type
       number             = thermostats_info%number_of_thermostats
       region             = simpar%part_region

       CALL thermostat_mapping_region(nhc%dis_type, natoms_local, nmol_local, const_mol,&
            tot_const, point, local_molecules, molecule_kind_set, simpar, error)

       ! Now we can allocate the target array s_kin and p_kin
       ALLOCATE (deg_of_freedom(number),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (nhc%s_kin(number),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (nhc%v_scale(number),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (nhc%p_kin(3,natoms_local),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (nhc%p_scale(3,natoms_local),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Allocate index array
       ALLOCATE ( nhc%index(number), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)

       ! Initialize deg_of_freedom
       deg_of_freedom = 0
       number         = 0
       IF ( have_global ) number = 1

       CALL massive_list_generate ( molecule_set, molecule_kind_set, &
            local_molecules, para_env, massive_atom_list, region, error )

       IF ( nhc % dis_type == do_thermo_no_communication ) THEN
          icount = 0
          IF ( region==do_region_global) THEN
             ! Global Region
             DO ikind = 1, nkind
                nmol_local =  local_molecules % n_el ( ikind )
                molecule_kind => molecule_kind_set ( ikind )
                DO imol_local = 1, nmol_local
                   icount = icount + 1
                   DO ii = 1, 3
                      DO jj = point ( 1, icount ), point ( 2, icount )
                         nhc%p_kin ( ii, jj ) % point => nhc % s_kin ( 1 )
                         nhc%p_scale ( ii, jj ) % point => nhc % v_scale ( 1 )
                      END DO
                   END DO
                   deg_of_freedom ( 1 ) = deg_of_freedom ( 1 ) + tot_const ( icount )
                   nhc%index(1) = 1
                END DO
             END DO
          ELSE IF ( region==do_region_molecule) THEN
             ! Molecular region
             DO ikind = 1, nkind
                nmol_local =  local_molecules % n_el ( ikind )
                molecule_kind => molecule_kind_set ( ikind )
                DO imol_local = 1, nmol_local
                   imol = local_molecules % list ( ikind ) % array ( imol_local )
                   icount = icount + 1
                   number = number + 1
                   nhc%index(number) = imol
                   deg_of_freedom ( number ) = const_mol ( icount )
                   DO jj = 1, 3
                      DO kk = point ( 1, icount ), point ( 2, icount )
                         nhc%p_kin(jj,kk) %point => nhc%s_kin(number)
                         nhc%p_scale(jj,kk) %point => nhc%v_scale(number)
                      END DO
                   END DO
                END DO
             END DO
          ELSE IF ( region==do_region_massive) THEN
             ! Massive Region
             DO ikind = 1, nkind
                nmol_local =  local_molecules % n_el ( ikind )
                molecule_kind => molecule_kind_set ( ikind )
                DO imol_local = 1, nmol_local
                   imol = local_molecules % list ( ikind ) % array ( imol_local )
                   molecule => molecule_set ( imol )
                   CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
                   icount = icount + 1
                   IF ( tot_const ( icount ) > 0 ) THEN 
                      CALL stop_program ( routineN, "Massive thermostats with constraints are impossible" )
                   END IF
                   atom_increment = 0
                   DO ii =  point ( 1, icount ), point ( 2, icount )
                      ipart = first_atom + atom_increment
                      DO  iii = 1, SIZE ( massive_atom_list )
                         IF ( ipart == massive_atom_list ( iii ) ) THEN
                            ielement = iii
                            EXIT
                         END IF
                      END DO
                      atom_increment = atom_increment + 1
                      DO jj = 1, 3
                         number = number + 1
                         nhc%index(number) = (ielement - 1)*3 + jj
                         nhc%p_kin(jj,ii) %point => nhc%s_kin(number)
                         nhc%p_scale(jj,ii) %point => nhc%v_scale(number)
                      END DO
                   END DO
                   IF ( first_atom + atom_increment -1 /= last_atom ) THEN
                      CALL stop_program ( routineN, "inconsistent mapping of particles" )
                   END IF
                END DO
             END DO
          ELSE
             CALL stop_program ( routineN, "Unknown region!" )
          END IF
       ELSE IF ( nhc % dis_type == do_thermo_communication ) THEN
          IF ( region==do_region_global) THEN
             ! Global Region
             DO ikind = 1, nkind
                molecule_kind => molecule_kind_set ( ikind )
                DO ii = 1, 3
                   DO jj = point ( 1, ikind ), point ( 2, ikind )
                      nhc%p_kin ( ii, jj ) % point => nhc % s_kin ( 1 )
                      nhc%p_scale ( ii, jj ) % point => nhc % v_scale ( 1 )
                   END DO
                END DO
                deg_of_freedom ( 1 ) = deg_of_freedom ( 1 ) + tot_const ( ikind )
                nhc%index(1) = 1
             END DO
          ELSE IF (region==do_region_molecule) THEN
             ! Molecular Region
             DO ikind = 1, nkind
                molecule_kind => molecule_kind_set ( ikind )
                CALL get_molecule_kind ( molecule_kind, natom = natom )
                number = number + 1
                nhc%index(number) = number
                deg_of_freedom ( number ) = deg_of_freedom ( number ) + tot_const ( ikind )
                DO jj = 1, 3
                   DO kk = point ( 1, ikind ), point ( 2, ikind )
                      nhc%p_kin ( jj, kk ) % point => nhc % s_kin ( number )
                      nhc%p_scale ( jj, kk ) % point => nhc % v_scale ( number )
                   END DO
                END DO
             END DO
          ELSE IF (region==do_region_massive) THEN
             ! Massive Region
             DO ikind = 1, nkind
                molecule_kind => molecule_kind_set ( ikind )
                IF ( tot_const ( ikind ) > 0 ) THEN
                   CALL stop_program ( routineN, "Massive thermostats with constraints are impossible" )
                END IF
                DO ii = point ( 1, ikind ), point ( 2, ikind )
                   DO jj = 1, 3
                      number = number + 1
                      nhc%index(number) = number
                      nhc%p_kin ( jj, ii ) % point => nhc % s_kin ( number )
                      nhc%p_scale ( jj, ii ) % point => nhc % v_scale ( number )
                   END DO
                END DO
             END DO
          ELSE
             CALL stop_program ( routineN, "Unknown region!" )
          ENDIF
       ELSE
          CALL stop_program ( routineN, "Unknown thermostat distribution option!" )
       END IF

       DEALLOCATE (const_mol,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (tot_const,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (point, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (nhc%nvt(nhc%nhc_len,number),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Sum up the number of degrees of freedom on each thermostat.
       ! first: initialize the target
       nhc%s_kin = 0.0_dp
       DO i = 1, 3
          DO j = 1, natoms_local
             nhc % p_kin ( i, j ) % point = nhc % p_kin ( i, j ) % point + 1
          END DO
       END DO

       ! if thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF ( nhc % dis_type == do_thermo_communication ) CALL mp_sum ( nhc % s_kin, para_env%group )

       !
       ! Now that we know how many there are stick this into nhc%nkt
       ! (number of degrees of freedom times k_B T for the first thermostat
       !  on the chain)
       nhc % num_nhc = number

       ! We know the total number of system thermostats.
       IF ( ( sum_of_thermostats == 1 ) .AND. ( nhc % dis_type /= do_thermo_no_communication ) ) THEN
          fac = nhc%s_kin(1)-deg_of_freedom(1)-simpar%nfree_rot_transl
          IF ( fac == 0.0_dp ) THEN 
             CALL stop_program ( routineN,'Zero degrees of freedom. Nothing to thermalize!' )
          END IF
          nhc % nvt ( 1, 1 ) % nkt = simpar%temp_ext*fac
       ELSE
          DO i = 1, number
             fac = ( nhc % s_kin ( i ) - deg_of_freedom ( i ) )
             nhc % nvt ( 1, i ) % nkt = simpar % temp_ext * fac 
          END DO
       END IF
       
       ! getting the number of degrees of freedom times k_B T for the rest
       ! of the chain
       DO i = 2, nhc % nhc_len
          nhc % nvt ( i, : ) % nkt = simpar % temp_ext
       END DO
       DEALLOCATE (deg_of_freedom,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( massive_atom_list, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE nhc_to_particle_mapping

!******************************************************************************
!!****** nhc_mapping/nhc_to_coef_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_coef_mapping
!!
!!   FUNCTION
!!     Creates the thermostatting for the
!!     polarization  coefficients'  degrees of freedom
!!
!!   AUTHOR
!!     CJM 10-Apr-01
!!
!!   MODIFICATION HISTORY
!!     gt 04-02-04 adapted to new kind of coefficents
!!     gt 11-02-05 adapted to new input
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE nhc_to_coef_mapping ( dyn_coeff_set, nhc, kg_coeff_section, para_env, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(section_vals_type), POINTER         :: kg_coeff_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_coef_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ikind, j, kg_coeff_ens, &
                                                kg_coeff_nhc, ncoef_local, &
                                                ncoef_of_kind, nkind, number, &
                                                stat
    LOGICAL                                  :: failure, nointer
    REAL(KIND=dp)                            :: kg_coeff_temp
    TYPE(section_vals_type), POINTER         :: nose_section

    failure = .FALSE.
    nose_section => section_vals_get_subs_vals(kg_coeff_section,"thermostat%nose",error=error)

    CALL section_vals_val_get(kg_coeff_section,"ensemble",i_val=kg_coeff_ens,error=error)
    CALL section_vals_val_get(kg_coeff_section,"temperature",r_val=kg_coeff_temp,error=error)

    SELECT CASE ( kg_coeff_ens )
    CASE  DEFAULT
       CALL stop_program(routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble )
       CALL section_vals_val_get(kg_coeff_section,"THERMOSTAT%REGION",i_val=kg_coeff_nhc,error=error)
       CALL section_vals_val_get(nose_section,"yoshida",i_val=nhc % nyosh ,error=error)
       CALL section_vals_val_get(nose_section,"mts",i_val=nhc % nc ,error=error)
       CALL section_vals_val_get(nose_section,"length",i_val=nhc % nhc_len ,error=error)
       number = 0
       ncoef_local = 0
       nointer = .TRUE.
       nkind = SIZE ( dyn_coeff_set%coeffs_of_kind)
       DO ikind=1,nkind
          IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
             CALL get_dyn_coeff(coeffs=dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                  ncoef_of_kind=ncoef_of_kind,error=error)
             ncoef_local= ncoef_local + ncoef_of_kind
          ENDIF
       END DO

       ! Counting the global number of thermostats
       IF ( kg_coeff_nhc == do_region_global )  THEN
          nointer=.FALSE.
          number=1
       ELSEIF ( kg_coeff_nhc == do_region_massive )  THEN
          number=ncoef_local
       ELSE
          CALL stop_program ( routineN, "Unknown region for KG COEFF!" )
       ENDIF

       ! Coefficients are local. communications needed only with GLOBAL
       IF ( nointer ) THEN
          nhc % dis_type = do_thermo_no_communication
       ELSE
          nhc % dis_type = do_thermo_communication
       END IF

       ! Now we can allocate the target array s_kin and p_kin
       ALLOCATE ( nhc %index(number), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % s_kin ( number ), STAT = stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % v_scale ( number  ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % p_kin ( 1, ncoef_local ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % p_scale ( 1, ncoef_local ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       nhc % num_nhc = number
       ! Begin the mapping loop
       i=0
       number = 0
       IF ( kg_coeff_nhc == do_region_global )  THEN
          DO i= 1, ncoef_local
             DO j = 1, SIZE ( nhc % p_kin, 1 )
                nhc % p_kin ( j, i ) % point => nhc % s_kin ( 1 )
                nhc % p_scale ( j, i ) % point => nhc % v_scale ( 1 )
             END DO
          END DO
          nhc%index(1) = 1
       ELSEIF (  kg_coeff_nhc == do_region_massive )  THEN
          DO i = 1, ncoef_local
             DO j = 1, SIZE ( nhc % p_kin, 1 )
                number = number + 1
                nhc % p_kin ( j, i ) % point => nhc % s_kin ( number )
                nhc % p_scale ( j, i ) % point => nhc % v_scale ( number )
                nhc%index(number) = number
             END DO
          END DO
       ENDIF

       ALLOCATE ( nhc % nvt ( nhc % nhc_len, nhc % num_nhc ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Find the number of degrees of freedom on each thermostat
       nhc % s_kin = 0.0_dp
       DO i = 1, ncoef_local
          DO j = 1, SIZE ( nhc % p_kin, 1 )
             nhc % p_kin ( j, i ) % point  = nhc % p_kin ( j, i ) % point  + 1.0_dp
          END DO
       ENDDO

       ! if thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF ( nhc % dis_type == do_thermo_communication ) CALL mp_sum ( nhc % s_kin, para_env%group )

       !
       ! Now that we know how many there are stick this into nhc % static % nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, nhc % num_nhc
          nhc % nvt ( 1, i ) % nkt = kg_coeff_temp * nhc % s_kin ( i )
       END DO

       ! getting the number of degrees of freedom times k_B T for the rest of the chain
       DO i = 2, nhc % nhc_len
          nhc % nvt ( i, : ) % nkt = kg_coeff_temp
       END DO
    END SELECT

  END SUBROUTINE nhc_to_coef_mapping

!******************************************************************************
!!****** nhc_mapping/nhc_to_shell_mapping [1.0] *
!!
!!   NAME
!!     nhc_to_shell_mapping
!!
!!   FUNCTION
!!     
!!
!!   AUTHOR
!!  
!!
!!   MODIFICATION HISTORY
!!   
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE nhc_to_shell_mapping(simpar, nhc, atomic_kind_set, local_particles, para_env,&
       error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_shell_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, icount, iii, ielement, iparticle_kind, iproc, isos, j, k, ncount,& 
      nkind, nparticle_local, nshell, nshell_local, number, offset, stat, sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: array_nshell_local
    INTEGER, DIMENSION(:), POINTER           :: massive_shell_list     
    INTEGER, DIMENSION(:), POINTER           :: part_index, work
    LOGICAL                                  :: failure, is_shell, nointer
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure = .FALSE.
    NULLIFY(array_nshell_local,atomic_kind,massive_shell_list,part_index,work)
    IF(simpar%shell_nose) THEN
       nhc%nyosh = simpar%shell_nyosh
       nhc%nc = simpar%shell_nc
       nhc%nhc_len = simpar%shell_nhclen 

       ! Counting the global number of thermostats
       sum_of_thermostats = 0
       ! variable to denote independent thermostats (no communication necessary)
       nointer = .TRUE.
       nshell_local = 0

       nkind =  SIZE(atomic_kind_set)
       DO iparticle_kind=1,nkind
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               shell_active=is_shell)
          IF(is_shell) THEN
             nparticle_local = local_particles%n_el(iparticle_kind)
             nshell_local = nshell_local  + nparticle_local
          END IF
       END DO


       IF(simpar%shell_region==do_region_global ) THEN
          ! Global Region
          nointer = .FALSE.
          sum_of_thermostats = 1
       ELSE IF (simpar%shell_region==do_region_molecule) THEN
          ! Molecule Regione
          sum_of_thermostats = nshell_local
       ELSE IF (simpar%shell_region==do_region_massive) THEN
          ! Massive Regione
          sum_of_thermostats = 3 * nshell_local
       ELSE
          CALL stop_program ( routineN, "Unknown region for SHELL!" )
       END IF
       number = sum_of_thermostats

       IF(nointer) THEN
          nhc % dis_type = do_thermo_no_communication
       ELSE
          nhc % dis_type = do_thermo_communication
       END IF

       ALLOCATE ( nhc %index(number), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % s_kin ( number ), STAT = stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % v_scale ( number  ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % p_kin ( 3,nshell_local),STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( nhc % p_scale ( 3,nshell_local),STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       nhc % num_nhc = number

       IF(nhc % dis_type == do_thermo_no_communication) THEN
         ! Particle mapping
         ALLOCATE (part_index(nshell_local),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         icount = 0
         DO iparticle_kind=1,nkind
            atomic_kind => atomic_kind_set(iparticle_kind)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                 shell_active=is_shell)
            nparticle_local = local_particles%n_el(iparticle_kind)
            IF(is_shell) THEN
               DO i = 1,nparticle_local
                 icount = icount + 1 
                 part_index (icount) = local_particles%list(iparticle_kind)%array(i)
               END DO 
            END IF
         END DO
         !!!!!!!
         ALLOCATE (array_nshell_local(para_env%num_pe), STAT=isos)
         CALL mp_allgather(nshell_local,array_nshell_local,para_env%group)
         nshell = SUM(array_nshell_local)
         ALLOCATE (massive_shell_list(nshell),STAT=isos)
         offset = 0
         DO iproc = 1,para_env%num_pe
           ncount = array_nshell_local(iproc)
           ALLOCATE (work(ncount), STAT = isos)
           IF(para_env%mepos == (iproc -1)) THEN
              DO i =1,ncount
                work(i) = part_index(i)
              END DO
           ELSE
              work(:) = 0
           END IF 
           CALL mp_bcast(work,iproc-1,para_env%group)
           DO i = 1,ncount
             massive_shell_list(offset+i) =work(i)
           END DO
           DEALLOCATE (work,STAT=isos)
           offset = offset + array_nshell_local(iproc)
         ENDDO
         DO i = 1,nshell
           DO j = i+1,nshell
             IF(massive_shell_list(i)>massive_shell_list(j)) THEN
                k=massive_shell_list(i)
                massive_shell_list(i) = massive_shell_list(j)
                massive_shell_list(j) = k
              END IF
            END DO
         END DO
         !!!!!!!
       END IF
       
       ! Begin the mapping loop
       i=0
       number = 0
       IF(simpar%shell_region==do_region_global ) THEN
          ! here distribute      
          DO i =1,nshell_local
             DO j = 1, 3
                nhc % p_kin ( j, i ) % point => nhc % s_kin ( 1 )
                nhc % p_scale ( j, i ) % point => nhc % v_scale ( 1 )
             END DO
          END DO
          nhc%index(1) = 1
       ELSEIF (simpar%shell_region==do_region_molecule) THEN
          DO i =1,nshell_local
             number = number + 1
             DO j = 1, 3
                nhc % p_kin ( j, i ) % point => nhc % s_kin ( number )
                nhc % p_scale ( j, i ) % point => nhc % v_scale ( number )
             END DO
             nhc % index ( number ) = part_index(i)
          END DO
          
       ELSEIF (simpar%shell_region==do_region_massive ) THEN 
          ! here distribute      
          DO i =1,nshell_local
             DO iii = 1, SIZE(massive_shell_list)
               IF ( part_index(i) == massive_shell_list(iii)) THEN
                  ielement = iii
                  EXIT
               END IF
             END DO
             DO j = 1, 3
                number = number + 1
                nhc % p_kin ( j, i ) % point => nhc % s_kin ( number )
                nhc % p_scale ( j, i ) % point => nhc % v_scale ( number )
                nhc % index ( number ) = (ielement -1) * 3 + j 
             END DO
          END DO
       END IF
       ALLOCATE ( nhc % nvt ( nhc % nhc_len, nhc % num_nhc ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Sum up the number of degrees of freedom on each thermostat.
       ! first: initialize the target
       nhc%s_kin = 0.0_dp
       DO j = 1, nshell_local
          DO i = 1, 3
             nhc % p_kin ( i, j ) % point = nhc % p_kin ( i, j ) % point + 1
          END DO
       END DO

       ! if thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF ( nhc % dis_type == do_thermo_communication ) CALL mp_sum ( nhc % s_kin, para_env%group )
       !
       ! Now that we know how many there are stick this into nhc%nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, nhc % num_nhc
          nhc % nvt ( 1, i ) % nkt = simpar%temp_sh_ext * nhc % s_kin ( i )
       END DO

       ! getting the number of degrees of freedom times k_B T for the rest of the chain
       DO i = 2, nhc % nhc_len
          nhc % nvt ( i, : ) % nkt = simpar%temp_sh_ext
       END DO

       IF(ASSOCIATED(part_index)) THEN
          DEALLOCATE(part_index, STAT=isos)
       END IF
       IF(ASSOCIATED(array_nshell_local)) THEN
          DEALLOCATE(array_nshell_local, STAT=isos)
       END IF
       IF(ASSOCIATED(massive_shell_list)) THEN
          DEALLOCATE(massive_shell_list, STAT=isos)
       END IF

    ELSE
       CALL stop_program(routineN, 'Never reach this point!' )
    END IF

  END SUBROUTINE nhc_to_shell_mapping

END MODULE nhc_mapping
