!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/thermostat_mapping [1.0] *
!!
!!   NAME
!!     thermostat_mapping
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - 10.2007 - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE thermostat_mapping

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: map_info_type
  USE input_constants,                 ONLY: do_region_global,&
                                             do_region_massive,&
                                             do_region_molecule,&
                                             do_thermo_communication,&
                                             do_thermo_no_communication
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PUBLIC :: thermostat_mapping_region,&
            init_baro_map_info,&
            init_coeff_map_info

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_mapping'

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** thermostat_mapping/thermostat_mapping_region [1.0] *
!!
!!   NAME
!!     thermostat_mapping_region
!!
!!   FUNCTION
!!     Main general setup thermostat regions (thermostat independent)
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE thermostat_mapping_region ( map_info, deg_of_freedom, massive_atom_list,&
       molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local, simpar, &
       number, region, gci, shell, error)

    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(OUT)                     :: natoms_local
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    INTEGER                                  :: number, region
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'thermostat_mapping_region', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nkind, nmol_local, &
                                                number_of_thermostats, stat
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    NULLIFY ( const_mol, tot_const, point)
    CPPostcondition(.NOT.ASSOCIATED(deg_of_freedom),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(massive_atom_list),cp_failure_level,routineP,error,failure)
    
    
    nkind = SIZE(molecule_kind_set)
    CALL mapping_region_evaluate(map_info%dis_type, natoms_local, nmol_local,&
         const_mol, tot_const, point, local_molecules, molecule_kind_set,&
         simpar, shell, error)

    ! Now we can allocate the target array s_kin and p_kin
    ALLOCATE (deg_of_freedom(number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%s_kin(number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%v_scale(number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%p_kin(3,natoms_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%p_scale(3,natoms_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Allocate index array
    ALLOCATE ( map_info%index(number), STAT = stat )
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    
    CALL massive_list_generate (molecule_set, molecule_kind_set, &
         local_molecules, para_env, massive_atom_list, region, shell, &
         error )
    
    CALL  thermostat_mapping_region_low(region, map_info, nkind, point,&
         deg_of_freedom, local_molecules, const_mol, massive_atom_list,&
         tot_const, molecule_set, number_of_thermostats, shell, gci,&
         error)
    
    check = (number==number_of_thermostats)
    CPPrecondition(check,cp_fatal_level,routineP,error,failure)
    DEALLOCATE (const_mol,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (tot_const,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (point, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE thermostat_mapping_region

!******************************************************************************
!!****** thermostat_mapping/thermostat_mapping_region_low [1.0] *
!!
!!   NAME
!!     thermostat_mapping_region_low
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
  SUBROUTINE  thermostat_mapping_region_low(region, map_info, nkind, point,&
       deg_of_freedom, local_molecules, const_mol, massive_atom_list, tot_const,&
       molecule_set, number, shell, gci, error)

    INTEGER, INTENT(IN)                      :: region
    TYPE(map_info_type), POINTER             :: map_info
    INTEGER                                  :: nkind
    INTEGER, DIMENSION(:, :), POINTER        :: point
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, DIMENSION(:), POINTER           :: const_mol, massive_atom_list, &
                                                tot_const
    TYPE(molecule_type), POINTER             :: molecule_set(:)
    INTEGER, INTENT(OUT)                     :: number
    LOGICAL, INTENT(IN)                      :: shell
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'thermostat_mapping_region_low', &
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
       check = ( map_info%dis_type == do_thermo_communication )
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       DO ikind = 1, nkind
          DO jj = point ( 1, ikind ), point ( 2, ikind )
             DO ii = 1, 3
                map_info%p_kin(ii,jj)%point => map_info%s_kin(1)
                map_info%p_scale(ii,jj)%point => map_info%v_scale(1)
             END DO
          END DO
          deg_of_freedom(1) = deg_of_freedom(1) + tot_const(ikind)
          map_info%index(1) = 1
          number = 1
       END DO
       deg_of_freedom(1) = deg_of_freedom(1) + nglob_cns
    ELSE IF ( region==do_region_molecule) THEN
       ! Molecular Region
       IF      ( map_info%dis_type == do_thermo_no_communication ) THEN
          ! This is the standard case.. 
          DO ikind = 1, nkind
             nmol_local =  local_molecules % n_el (ikind)
             DO imol_local = 1, nmol_local
                imol = local_molecules%list(ikind)%array(imol_local)
                icount = icount + 1
                number = number + 1
                map_info%index(number) = imol
                deg_of_freedom ( number ) = const_mol ( icount )
                DO kk = point ( 1, icount ), point ( 2, icount )
                   DO jj = 1, 3
                      map_info%p_kin(jj,kk) %point => map_info%s_kin(number)
                      map_info%p_scale(jj,kk) %point => map_info%v_scale(number)
                   END DO
                END DO
             END DO
          END DO
       ELSE IF ( map_info%dis_type == do_thermo_communication ) THEN
          ! This case is quite rare and happens only when we have one molecular
          ! kind and one molecule..
          DO ikind = 1, nkind
             number = number + 1
             map_info%index(number) = number
             deg_of_freedom ( number ) = deg_of_freedom ( number ) + tot_const ( ikind )
             DO kk = point ( 1, ikind ), point ( 2, ikind )
                DO jj = 1, 3
                   map_info%p_kin ( jj, kk ) % point => map_info%s_kin ( number )
                   map_info%p_scale ( jj, kk ) % point => map_info%v_scale ( number )
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
       check = ( map_info%dis_type == do_thermo_no_communication )
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       DO ikind = 1, nkind
          nmol_local =  local_molecules % n_el ( ikind )
          DO imol_local = 1, nmol_local
             icount = icount + 1
             imol = local_molecules % list ( ikind ) % array ( imol_local )
             molecule => molecule_set ( imol )
             CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom,&
                  first_shell=first_shell, last_shell=last_shell)
             IF (shell) THEN
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
                   map_info%index(number) = (ielement - 1)*3 + jj
                   map_info%p_kin(jj,ii) %point => map_info%s_kin(number)
                   map_info%p_scale(jj,ii) %point => map_info%v_scale(number)
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

  END SUBROUTINE thermostat_mapping_region_low

!******************************************************************************
!!****** thermostat_mapping/mapping_region_evaluate [1.0] *
!!
!!   NAME
!!     mapping_region_evaluate
!!
!!   FUNCTION
!!     creates the mapping between the system and the thermostats
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - 10.2007 - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE mapping_region_evaluate(dis_type, natoms_local, nmol_local, const_mol,&
       tot_const, point, local_molecules, molecule_kind_set, simpar, shell,&
       error)
    INTEGER, INTENT(IN)                      :: dis_type
    INTEGER, INTENT(OUT)                     :: natoms_local, nmol_local
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    LOGICAL, INTENT(IN)                      :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mapping_region_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: atm_offset, icount, ikind, ilist, imol, katom, natom, nc, &
      nfixd, nfixd_restraint, nkind, nmol_per_kind, nmolecule, nshell, stat
    LOGICAL                                  :: failure
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    natoms_local = 0
    nmol_local   = 0
    nkind = SIZE ( molecule_kind_set )
    ! Compute the TOTAL number of molecules and atoms on THIS PROC and
    ! TOTAL number of molecules of IKIND on THIS PROC
    DO ikind = 1, nkind
       molecule_kind => molecule_kind_set ( ikind )
       CALL get_molecule_kind ( molecule_kind, natom=natom, nshell=nshell )
       IF (shell) THEN
          IF (nshell/=0) THEN
             natoms_local = natoms_local + nshell * local_molecules % n_el ( ikind )
             nmol_local   = nmol_local + local_molecules % n_el ( ikind )
          END IF
       ELSE
          natoms_local = natoms_local + natom * local_molecules % n_el ( ikind )
          nmol_local   = nmol_local + local_molecules % n_el ( ikind )
       END IF
    END DO

    CPPostcondition(.NOT.ASSOCIATED(const_mol),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(tot_const),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(point),cp_failure_level,routineP,error,failure)
    IF ( dis_type == do_thermo_no_communication ) THEN
       ALLOCATE ( const_mol ( nmol_local ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ALLOCATE ( tot_const ( nmol_local ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ALLOCATE ( point ( 2, nmol_local ), STAT = stat )

       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       point (:,:)= 0
       atm_offset = 0
       icount = 0
       DO ikind = 1, nkind
          nmol_per_kind =  local_molecules % n_el ( ikind )
          molecule_kind => molecule_kind_set ( ikind )
          CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
               fixd_list=fixd_list, nshell=nshell)
          IF (shell) natom = nshell
          DO imol = 1, nmol_per_kind
             icount = icount + 1
             point ( 1, icount )  = atm_offset + 1
             point ( 2, icount )  = atm_offset + natom
             IF (.NOT.shell) THEN
                ! nc keeps track of all constraints but not fixed ones..
                ! Let's identify fixed atoms for this molecule
                nfixd = 0
                IF (ASSOCIATED(fixd_list)) THEN
                   DO katom = point ( 1, icount ), point ( 2, icount )
                      DO ilist = 1, SIZE(fixd_list)
                         IF ( ( katom == fixd_list(ilist)%fixd        ) .AND. &
                              (.NOT. fixd_list(ilist)%restraint%active)) nfixd=nfixd+1
                      END DO
                   END DO
                END IF
                const_mol ( icount ) = nc + 3*nfixd
                tot_const ( icount ) = const_mol ( icount )
             END IF
             atm_offset =  point ( 2, icount )
          END DO
       END DO
    ELSE IF ( dis_type == do_thermo_communication ) THEN
       ALLOCATE ( const_mol ( nkind ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( tot_const ( nkind ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( point ( 2, nkind ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       point ( :, : ) = 0
       atm_offset = 0
       ! nc keeps track of all constraints but not fixed ones..
       DO ikind = 1, nkind
          nmol_per_kind =  local_molecules % n_el ( ikind )
          molecule_kind => molecule_kind_set ( ikind )
          CALL get_molecule_kind ( molecule_kind, nconstraint=nc, natom = natom,&
               nmolecule=nmolecule, nfixd=nfixd, nfixd_restraint=nfixd_restraint,&
               nshell=nshell)
          IF (shell) natom = nshell
          IF (.NOT.shell) THEN
             const_mol ( ikind ) = nc
             ! Let's consider the fixed atoms only for the total number of constraints
             ! in case we are in REPLICATED/INTERACTING thermostats
             tot_const ( ikind ) = const_mol ( ikind ) * nmolecule + 3*(nfixd-nfixd_restraint)
          END IF
          point ( 1, ikind ) = atm_offset + 1
          point ( 2, ikind ) = atm_offset + natom * nmol_per_kind
          atm_offset =  point ( 2, ikind )
       END DO
    ENDIF
    IF (( .NOT. simpar % constraint ).OR.shell) THEN
       const_mol = 0.0_dp
       tot_const = 0.0_dp
    END IF

  END SUBROUTINE mapping_region_evaluate

!******************************************************************************
!!****** thermostat_mapping/massive_list_generate [1.0] *
!!
!!   NAME
!!    massive_list_generate 
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
  SUBROUTINE massive_list_generate ( molecule_set, molecule_kind_set, &
       local_molecules, para_env, massive_atom_list, region, shell, error )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, POINTER                         :: massive_atom_list( : )
    INTEGER, INTENT(IN)                      :: region
    LOGICAL, INTENT(IN)                      :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'massive_list_generate', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, first_shell, i, ikind, imol, iproc, j, k, natom, &
      ncount, nkind, nmol_per_kind, nshell, num_massive_atm, &
      num_massive_atm_local, offset, stat
    INTEGER, DIMENSION(:), POINTER           :: array_num_massive_atm, &
                                                local_atm_list, work
    LOGICAL                                  :: failure
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    num_massive_atm_local = 0
    NULLIFY(local_atm_list)
    CALL reallocate(local_atm_list,1,num_massive_atm_local)

    nkind = SIZE ( molecule_kind_set )
    DO ikind = 1,  nkind
       nmol_per_kind =  local_molecules%n_el(ikind)
       DO imol = 1, nmol_per_kind
          i = local_molecules%list(ikind)%array(imol)
          molecule => molecule_set ( i )
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind(molecule_kind,natom=natom,nshell=nshell)
          IF(region == do_region_massive) THEN
             IF (shell) THEN
                natom = nshell
             END IF
             num_massive_atm_local = num_massive_atm_local + natom
             CALL reallocate(local_atm_list,1,num_massive_atm_local)
             CALL get_molecule (molecule,first_atom=first_atom,first_shell=first_shell)
             IF (shell) THEN
                first_atom = first_shell
             END IF
             DO j=1,natom
                local_atm_list(num_massive_atm_local-natom+j) = first_atom -1 + j
             END DO
          END IF
       END DO
    END DO

    ALLOCATE(array_num_massive_atm(para_env%num_pe), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL mp_allgather(num_massive_atm_local,array_num_massive_atm,para_env%group)

    num_massive_atm = SUM(array_num_massive_atm)
    ALLOCATE(massive_atom_list(num_massive_atm), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    offset = 0
    DO iproc=1,para_env%num_pe
       ncount = array_num_massive_atm(iproc)
       ALLOCATE(work(ncount), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF(para_env%mepos == (iproc-1)) THEN
          DO i=1,ncount
             work(i) = local_atm_list(i)
          END DO
       ELSE
          work(:) = 0
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       DO i=1,ncount
          massive_atom_list(offset+i) = work(i)
       END DO
       DEALLOCATE(work, STAT = stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       offset = offset + array_num_massive_atm(iproc)
    END DO

    DO i=1,num_massive_atm
       DO j=1,num_massive_atm
          IF(massive_atom_list(i)<massive_atom_list(j)) THEN
             k = massive_atom_list ( i )
             massive_atom_list ( i )  = massive_atom_list ( j )
             massive_atom_list ( j ) = k
          ENDIF
       END DO
    END DO

    DEALLOCATE(local_atm_list, STAT = stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(array_num_massive_atm, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE massive_list_generate

!******************************************************************************
!!****** thermostat_mapping/init_baro_map_info [1.0] *
!!
!!   NAME
!!     init_baro_map_info
!!
!!   FUNCTION
!!     Initialize the map_info for barostat thermostat
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - 10.2007 - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE  init_baro_map_info(map_info, number, num_thermo, error)
    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, INTENT(IN)                      :: number, num_thermo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_baro_map_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    ALLOCATE (map_info%s_kin(num_thermo),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%v_scale(num_thermo),stat=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%p_kin(1,number),stat=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (map_info%p_scale(1,number),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Allocate the index array
    ALLOCATE (map_info%index(1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    ! Begin the mapping loop
    DO i = 1, number
       map_info%p_kin(1,i)%point   => map_info%s_kin(1)
       map_info%p_scale(1,i)%point => map_info%v_scale(1)
    END DO
    map_info%index(1) = 1

  END SUBROUTINE init_baro_map_info


!******************************************************************************
!!****** thermostat_mapping/init_coeff_map_info [1.0] *
!!
!!   NAME
!!     init_coeff_map_info
!!
!!   FUNCTION
!!     Initialize the map_info for coefficients thermostat
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - 10.2007 - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE  init_coeff_map_info(map_info, number, kg_coeff_nhc, ncoef_local, error)
    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, INTENT(IN)                      :: number, kg_coeff_nhc, &
                                                ncoef_local
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_coeff_map_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, loc_number, stat
    LOGICAL                                  :: failure

    ALLOCATE(map_info%index(number),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_info%s_kin (number),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_info%v_scale(number),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_info%p_kin(1,ncoef_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(map_info%p_scale(1,ncoef_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Begin the mapping loop
    i=0
    loc_number = 0
    IF ( kg_coeff_nhc == do_region_global )  THEN
       DO i= 1, ncoef_local
          DO j = 1, SIZE(map_info%p_kin, 1)
             map_info%p_kin(j,i)%point   => map_info%s_kin(1)
             map_info%p_scale(j,i)%point => map_info%v_scale(1)
          END DO
       END DO
       map_info%index(1) = 1
    ELSEIF (  kg_coeff_nhc == do_region_massive )  THEN
       DO i = 1, ncoef_local
          DO j = 1, SIZE(map_info%p_kin, 1)
             loc_number = loc_number + 1
             map_info%p_kin(j,i)%point   => map_info%s_kin(loc_number)
             map_info%p_scale(j,i)%point => map_info%v_scale(loc_number)
             map_info%index(loc_number) = loc_number
          END DO
       END DO
    ENDIF
  END SUBROUTINE init_coeff_map_info

END MODULE thermostat_mapping

!******************************************************************************
