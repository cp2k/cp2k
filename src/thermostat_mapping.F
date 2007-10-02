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
  USE input_constants,                 ONLY: do_region_massive,&
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
                                             molecule_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PUBLIC :: thermostat_mapping_region, massive_list_generate

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_mapping'

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** thermostat_mapping/thermo_mapping_region [1.0] *
!!
!!   NAME
!!     thermo_mapping_region
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
  SUBROUTINE thermostat_mapping_region(dis_type, natoms_local, nmol_local, const_mol,&
       tot_const, point, local_molecules, molecule_kind_set, simpar, error)
    INTEGER, INTENT(IN)                      :: dis_type
    INTEGER, INTENT(OUT)                     :: natoms_local, nmol_local
    INTEGER, DIMENSION(:), POINTER           :: const_mol, tot_const
    INTEGER, DIMENSION(:, :), POINTER        :: point
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'thermostat_mapping_region', &
      routineP = moduleN//':'//routineN

    INTEGER :: atm_offset, icount, ikind, ilist, imol, katom, natom, nc, &
      nfixd, nfixd_restraint, nkind, nmol_per_kind, nmolecule, stat
    LOGICAL                                  :: failure
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    natoms_local = 0
    nmol_local = 0
    nkind = SIZE ( molecule_kind_set )
    ! Compute the TOTAL number of molecules and atoms on THIS PROC and
    ! TOTAL number of molecules of IKIND on THIS PROC
    DO ikind = 1, nkind
       molecule_kind => molecule_kind_set ( ikind )
       CALL get_molecule_kind ( molecule_kind, natom=natom )
       natoms_local = natoms_local + natom * local_molecules % n_el ( ikind )
       nmol_local   = nmol_local + local_molecules % n_el ( ikind )
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
               fixd_list=fixd_list)
          ! nc keeps track of all constraints but not fixed ones..
          DO imol = 1, nmol_per_kind
             icount = icount + 1
             point ( 1, icount )  = atm_offset + 1
             point ( 2, icount )  = atm_offset + natom
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
               nmolecule=nmolecule, nfixd=nfixd,&
               nfixd_restraint=nfixd_restraint)
          const_mol ( ikind ) = nc
          ! Let's consider the fixed atoms only for the total number of constraints
          ! in case we are in REPLICATED/INTERACTING thermostats
          tot_const ( ikind ) = const_mol ( ikind ) * nmolecule + 3*(nfixd-nfixd_restraint)
          point ( 1, ikind ) = atm_offset + 1
          point ( 2, ikind ) = atm_offset + natom * nmol_per_kind
          atm_offset =  point ( 2, ikind )
       END DO
    ENDIF
    IF ( .NOT. simpar % constraint ) THEN
       const_mol = 0.0_dp
       tot_const = 0.0_dp
    END IF

  END SUBROUTINE thermostat_mapping_region

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
       local_molecules, para_env, massive_atom_list, region, error )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, POINTER                         :: massive_atom_list( : )
    INTEGER, INTENT(IN)                      :: region
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'massive_list_generate', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, i, ikind, imol, iproc, j, k, natom, ncount, nkind, &
      nmol_per_kind, num_massive_atm, num_massive_atm_local, offset, stat
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
          CALL get_molecule_kind(molecule_kind,natom=natom)
          IF(region == do_region_massive) THEN
             num_massive_atm_local = num_massive_atm_local + natom
             CALL reallocate(local_atm_list,1,num_massive_atm_local)
             CALL get_molecule (molecule,first_atom=first_atom)
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

END MODULE thermostat_mapping

!******************************************************************************
