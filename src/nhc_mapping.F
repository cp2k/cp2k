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

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             get_dyn_coeff
  USE extended_system_types,           ONLY: debug_isotropic_limit,&
                                             lnhc_parameters_type
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
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE termination,                     ONLY: stop_program
  USE thermostat_mapping,              ONLY: massive_list_generate,&
                                             thermostat_mapping_region
  USE thermostat_types,                ONLY: thermostat_info_type
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
  SUBROUTINE nhc_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, nhc, para_env, gci, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_particle_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, natoms_local, stat, &
                                                sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fac

    failure = .FALSE.
    NULLIFY (massive_atom_list, deg_of_freedom)
    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program ( routineN, 'Unknown ensemble!' )
    CASE ( nve_ensemble, isokin_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program ( routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble )
       
       CALL setup_nhc_thermostat(nhc, thermostat_info, simpar%nyosh, simpar%nc,&
            simpar%nhclen, simpar%part_region, deg_of_freedom, massive_atom_list,&
            molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, &
            simpar, sum_of_thermostats, gci, error=error)
       
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

       ! We know the total number of system thermostats.
       IF ( ( sum_of_thermostats == 1 ) .AND. ( nhc % dis_type /= do_thermo_no_communication ) ) THEN
          fac = nhc%s_kin(1)-deg_of_freedom(1)-simpar%nfree_rot_transl
          IF ( fac == 0.0_dp ) THEN 
             CALL stop_program ( routineN,'Zero degrees of freedom. Nothing to thermalize!' )
          END IF
          nhc % nvt ( 1, 1 ) % nkt = simpar%temp_ext*fac
       ELSE
          DO i = 1, nhc % num_nhc
             fac = ( nhc % s_kin ( i ) - deg_of_freedom ( i ) )
             nhc % nvt ( 1, i ) % nkt = simpar % temp_ext * fac 
          END DO
       END IF
       
       ! Getting the number of degrees of freedom times k_B T for the rest
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
!!****** nhc_mapping/setup_nhc_thermostat [1.0] *
!!
!!   NAME
!!     setup_nhc_thermostat
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
  SUBROUTINE setup_nhc_thermostat ( nhc, thermostat_info, nyosh, nc, nhclen,&
       part_region, deg_of_freedom, massive_atom_list, molecule_kind_set,&
       local_molecules, molecule_set, para_env, natoms_local, simpar, sum_of_thermostats,&
       gci, shell, error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    INTEGER, INTENT(IN)                      :: nyosh, nc, nhclen, part_region
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_nhc_thermostat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nkind, nmol_local, number, &
                                                number_of_thermostats, &
                                                region, stat
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    LOGICAL                                  :: check, do_shell, failure

    failure = .FALSE.
    NULLIFY ( const_mol, tot_const, point)
    CPPostcondition(.NOT.ASSOCIATED(deg_of_freedom),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(massive_atom_list),cp_failure_level,routineP,error,failure)
    do_shell = .FALSE.
    IF (PRESENT(shell)) do_shell = shell

    nhc%nyosh = nyosh
    nhc%nc = nc
    nhc%nhc_len = nhclen
    
    nkind = SIZE(molecule_kind_set)
    sum_of_thermostats = thermostat_info%sum_of_thermostats
    nhc%dis_type       = thermostat_info%dis_type
    number             = thermostat_info%number_of_thermostats
    region             = part_region

    CALL thermostat_mapping_region(nhc%dis_type, natoms_local, nmol_local,&
         const_mol, tot_const, point, local_molecules, molecule_kind_set,&
         simpar, do_shell, error)

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
    
    CALL massive_list_generate (molecule_set, molecule_kind_set, &
         local_molecules, para_env, massive_atom_list, region, do_shell, &
         error )
    
    CALL nhc_region_mapping (region, nhc, nkind, point, deg_of_freedom,&
         local_molecules, const_mol, massive_atom_list, tot_const,&
         molecule_set, number_of_thermostats, do_shell, gci, error)
    
    check = (number==number_of_thermostats)
    CPPrecondition(check,cp_fatal_level,routineP,error,failure)
    DEALLOCATE (const_mol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (tot_const,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (point, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (nhc%nvt(nhc%nhc_len,number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Now that we know how many there are stick this into nhc%nkt
    ! (number of degrees of freedom times k_B T for the first thermostat
    !  on the chain)
    nhc % num_nhc = number

  END SUBROUTINE setup_nhc_thermostat

!******************************************************************************
!!****** nhc_mapping/nhc_region_mapping [1.0] *
!!
!!   NAME
!!     nhc_region_mapping
!!
!!   FUNCTION
!!     Performs the real mapping for the thermostat region
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE nhc_region_mapping (region, nhc, nkind, point, deg_of_freedom, local_molecules,&
       const_mol, massive_atom_list, tot_const, molecule_set, number, do_shell, gci, error)

    INTEGER, INTENT(IN)                      :: region
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    INTEGER                                  :: nkind
    INTEGER, DIMENSION(:, :), POINTER        :: point
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, DIMENSION(:), POINTER           :: const_mol, massive_atom_list, &
                                                tot_const
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    INTEGER, INTENT(OUT)                     :: number
    LOGICAL, INTENT(IN)                      :: do_shell
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_region_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, first_shell, icount, ielement, ii, iii, ikind, &
      imol, imol_local, ipart, jj, k, kk, last_atom, last_shell, nglob_cns, &
      nmol_local
    LOGICAL                                  :: check, failure, &
                                                global_constraints
    TYPE(molecule_type), POINTER             :: molecule

    failure        = .FALSE.
    global_constraints = ASSOCIATED(gci)
    deg_of_freedom = 0
    icount         = 0
    number         = 0
    nglob_cns      = 0
    IF (global_constraints) nglob_cns = gci%ntot-gci%nrestraint
    IF ( region==do_region_global) THEN
       ! Global Region
       check = ( nhc % dis_type == do_thermo_communication )
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       DO ikind = 1, nkind
          DO jj = point ( 1, ikind ), point ( 2, ikind )
             DO ii = 1, 3
                nhc%p_kin(ii,jj)%point => nhc%s_kin(1)
                nhc%p_scale(ii,jj)%point => nhc%v_scale(1)
             END DO
          END DO
          deg_of_freedom(1) = deg_of_freedom(1) + tot_const(ikind)
          nhc%index(1) = 1
          number = 1
       END DO
       deg_of_freedom(1) = deg_of_freedom(1) + nglob_cns
    ELSE IF ( region==do_region_molecule) THEN
       ! Molecular Region
       IF      ( nhc % dis_type == do_thermo_no_communication ) THEN
          ! This is the standard case.. 
          DO ikind = 1, nkind
             nmol_local =  local_molecules % n_el (ikind)
             DO imol_local = 1, nmol_local
                imol = local_molecules%list(ikind)%array(imol_local)
                icount = icount + 1
                number = number + 1
                nhc%index(number) = imol
                deg_of_freedom ( number ) = const_mol ( icount )
                DO kk = point ( 1, icount ), point ( 2, icount )
                   DO jj = 1, 3
                      nhc%p_kin(jj,kk) %point => nhc%s_kin(number)
                      nhc%p_scale(jj,kk) %point => nhc%v_scale(number)
                   END DO
                END DO
             END DO
          END DO
       ELSE IF ( nhc % dis_type == do_thermo_communication ) THEN
          ! This case is quite rare and happens only when we have one molecular
          ! kind and one molecule..
          DO ikind = 1, nkind
             number = number + 1
             nhc%index(number) = number
             deg_of_freedom ( number ) = deg_of_freedom ( number ) + tot_const ( ikind )
             DO kk = point ( 1, ikind ), point ( 2, ikind )
                DO jj = 1, 3
                   nhc%p_kin ( jj, kk ) % point => nhc % s_kin ( number )
                   nhc%p_scale ( jj, kk ) % point => nhc % v_scale ( number )
                END DO
             END DO
          END DO
       ELSE
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       IF (nglob_cns/=0) THEN
          CALL stop_program ( routineN, "Molecular thermostats with global constraints are impossible!" )
       END IF
    ELSE IF ( region==do_region_massive) THEN
       ! Massive Region
       check = ( nhc % dis_type == do_thermo_no_communication )
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       DO ikind = 1, nkind
          nmol_local =  local_molecules % n_el ( ikind )
          DO imol_local = 1, nmol_local
             icount = icount + 1
             imol = local_molecules % list ( ikind ) % array ( imol_local )
             molecule => molecule_set ( imol )
             CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom,&
                  first_shell=first_shell, last_shell=last_shell)
             IF (do_shell) THEN
                first_atom = first_shell
                last_atom  = last_shell
             ELSE
                IF ((tot_const(icount)>0).OR.(nglob_cns/=0)) THEN 
                   CALL stop_program ( routineN, "Massive thermostats with constraints are impossible!" )
                END IF
             END IF
             k = 0
             DO ii =  point ( 1, icount ), point ( 2, icount )
                ipart = first_atom + k
                DO  iii = 1, SIZE ( massive_atom_list )
                   IF ( ipart == massive_atom_list ( iii ) ) THEN
                      ielement = iii
                      EXIT
                   END IF
                END DO
                k = k + 1
                DO jj = 1, 3
                   number = number + 1
                   nhc%index(number) = (ielement - 1)*3 + jj
                   nhc%p_kin(jj,ii) %point => nhc%s_kin(number)
                   nhc%p_scale(jj,ii) %point => nhc%v_scale(number)
                END DO
             END DO
             IF ( first_atom + k -1 /= last_atom ) THEN
                CALL stop_program ( routineN, "inconsistent mapping of particles" )
             END IF
          END DO
       END DO
    ELSE
       CALL stop_program ( routineN, "Unknown region!" )
    END IF

  END SUBROUTINE nhc_region_mapping

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
  SUBROUTINE nhc_to_shell_mapping( thermostat_info, simpar, local_molecules,&
       molecule_set, molecule_kind_set, nhc, para_env, gci, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nhc_to_shell_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nshell_local, stat, &
                                                sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_shell_list
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY (massive_shell_list, deg_of_freedom)

    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program ( routineN, 'Unknown ensemble!' )
    CASE ( nve_ensemble, isokin_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program ( routineN, 'Never reach this point!' )
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble )

       CALL setup_nhc_thermostat(nhc, thermostat_info, simpar%shell_nyosh, simpar%shell_nc,&
            simpar%shell_nhclen, simpar%shell_region, deg_of_freedom, massive_shell_list,&
            molecule_kind_set, local_molecules, molecule_set, para_env, nshell_local, &
            simpar, sum_of_thermostats, gci, shell=.TRUE., error=error)

       ! Sum up the number of degrees of freedom on each thermostat.
       ! first: initialize the target
       nhc%s_kin = 0.0_dp
       DO j = 1, nshell_local
          DO i = 1, 3
             nhc % p_kin ( i, j ) % point = nhc % p_kin ( i, j ) % point + 1
          END DO
       END DO

       ! If thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF ( nhc % dis_type == do_thermo_communication ) CALL mp_sum ( nhc % s_kin, para_env%group )

       ! Now that we know how many there are stick this into nhc%nkt
       ! (number of degrees of freedom times k_B T )
       DO i = 1, nhc % num_nhc
          nhc % nvt ( 1, i ) % nkt = simpar%temp_sh_ext * nhc % s_kin ( i )
       END DO

       ! Getting the number of degrees of freedom times k_B T for the rest of the chain
       DO i = 2, nhc % nhc_len
          nhc % nvt ( i, : ) % nkt = simpar%temp_sh_ext
       END DO
       DEALLOCATE (deg_of_freedom,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE ( massive_shell_list, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE nhc_to_shell_mapping

END MODULE nhc_mapping
