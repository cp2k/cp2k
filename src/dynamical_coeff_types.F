!---------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/dynamical_coeff_types [1.0] *
!!
!!   NAME
!!     dynamical_coeff_type
!!   FUNCTION
!!     Define the data structure for the dinamical coefficents .
!!
!!   AUTHOR
!!     gloria (09.12.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE dynamical_coeff_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dyn_coeff_mass_atom,&
                                             dyn_coeff_mass_coeff,&
                                             dyn_coeff_mass_global,&
                                             dyn_coeff_mass_kind,&
                                             kg_coeff_init_start,&
                                             kg_coeff_pos_restart,&
                                             kg_coeff_pos_vel_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE termination,                     ONLY: stop_memory
  USE util,                            ONLY: gasdev,&
                                             locate

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='dynamical_coeff_type'
  INTEGER, PRIVATE, SAVE :: last_coeff_id=0
  INTEGER, PRIVATE, SAVE :: last_coeff_set_id=0

  INTEGER, PARAMETER,PUBLIC :: dyn_coeff_replicated=1, dyn_coeff_distributed=2

  TYPE dyn_coeff_set_type
   INTEGER :: ref_count, id_nr
   INTEGER :: mass_storage, distribution_method
   INTEGER :: ncoef_atom_max, ncoef_tot, nel_tot
   LOGICAL :: propagate
   REAL(KIND=dp) :: global_mass
   TYPE(distribution_1d_type), POINTER :: distribution
   TYPE(dyn_coeff_p_type), POINTER, DIMENSION(:) :: coeffs_of_kind
  END TYPE dyn_coeff_set_type

  TYPE dyn_coeff_set_p_type
   TYPE(dyn_coeff_set_type), POINTER :: dyn_coeff_set
  END TYPE dyn_coeff_set_p_type

  TYPE dyn_coeff_type
   INTEGER :: ref_count, id_nr
   INTEGER :: n_els
   INTEGER :: ncoef_atom
   REAL(KIND=dp), DIMENSION(:,:), POINTER :: pos, masses, &
     vel, forces
  END TYPE dyn_coeff_type

  TYPE dyn_coeff_p_type
   TYPE(dyn_coeff_type), POINTER :: coeffs
  END TYPE dyn_coeff_p_type

! *** Public structures ***

  PUBLIC :: dyn_coeff_set_type, dyn_coeff_set_p_type,&
            dyn_coeff_type, dyn_coeff_p_type

! *** Public Subroutines ***
 
  PUBLIC :: dyn_coeff_set_create, dyn_coeff_set_dump, dyn_coeff_set_initialize,&
            dyn_coeff_set_release, dyn_coeff_set_retain,&
            get_dyn_coeff_set, dyn_coeff_create, setup_dyn_coeff, &
            dyn_coeff_release, get_dyn_coeff, read_coeff_restart,&
            write_coeff_restart


! *****************************************************************************
  CONTAINS

! *****************************************************************************

!!*** **********************************************************************
SUBROUTINE dyn_coeff_set_create(dyn_coeff_set,atomic_kind_set,&
                                distribution, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      pointer                                :: atomic_kind_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: distribution
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, n_els, nkind, stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set

  failure=.FALSE.

  NULLIFY(atomic_kind,aux_basis_set)
  ALLOCATE(dyn_coeff_set, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_coeff_set_id=last_coeff_set_id+1
     dyn_coeff_set%id_nr=last_coeff_set_id
     dyn_coeff_set%ref_count=1
     NULLIFY(dyn_coeff_set%distribution, dyn_coeff_set%coeffs_of_kind)
!    set coeffs_of_kind
     nkind=size(atomic_kind_set)
     ALLOCATE(dyn_coeff_set%coeffs_of_kind(nkind), stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ikind=1,nkind
       NULLIFY(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind,aux_basis_set=aux_basis_set)
       IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
       CALL dyn_coeff_create(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                             error=error)
     END DO
     IF (PRESENT(distribution)) THEN
       IF (ASSOCIATED(distribution)) THEN
         dyn_coeff_set%distribution_method=dyn_coeff_distributed
         DO ikind=1,nkind 
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind,aux_basis_set=aux_basis_set)
            IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
            n_els = distribution%n_el(ikind)
            CALL setup_dyn_coeff(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
               aux_basis_set=aux_basis_set, n_els=n_els, error=error)
         END DO
         CALL distribution_1d_retain(distribution,error=error)
         CALL distribution_1d_release(dyn_coeff_set%distribution,error=error)
         dyn_coeff_set%distribution => distribution
       END IF
     ELSE 
       dyn_coeff_set%distribution_method=dyn_coeff_replicated
       DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              aux_basis_set=aux_basis_set,&
                              natom=n_els)         
         IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
         CALL setup_dyn_coeff(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
               aux_basis_set=aux_basis_set, n_els=n_els, error=error)
       END DO
     END IF
  END IF

END SUBROUTINE dyn_coeff_set_create

!***************************************************************************
SUBROUTINE dyn_coeff_set_retain(dyn_coeff_set, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(dyn_coeff_set),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(dyn_coeff_set%ref_count>0,cp_failure_level,routineP,error)
     dyn_coeff_set%ref_count=dyn_coeff_set%ref_count+1
  END IF

END SUBROUTINE dyn_coeff_set_retain

!!*** **********************************************************************

SUBROUTINE dyn_coeff_set_release(dyn_coeff_set, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(dyn_coeff_set)) THEN
     CPPreconditionNoFail(dyn_coeff_set%ref_count>0,cp_failure_level,routineP,error)
     dyn_coeff_set%ref_count=dyn_coeff_set%ref_count-1
     IF (dyn_coeff_set%ref_count==0) THEN
        CALL distribution_1d_release(dyn_coeff_set%distribution, error=error)
        IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind)) THEN
           DO i=1,SIZE(dyn_coeff_set%coeffs_of_kind)
              IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(i)%coeffs)) THEN
                 CALL dyn_coeff_release(dyn_coeff_set%coeffs_of_kind(i)%coeffs,error=error)
              ENDIF
           END DO
           DEALLOCATE(dyn_coeff_set%coeffs_of_kind,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(dyn_coeff_set, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(dyn_coeff_set)

END SUBROUTINE dyn_coeff_set_release

!***************************************************************************

SUBROUTINE get_dyn_coeff_set(dyn_coeff_set,distribution, coeffs_of_kind,&
                             ncoefs,mass_storage,global_mass,id_nr,ref_count,error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: distribution
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: coeffs_of_kind
    INTEGER, INTENT(out), OPTIONAL           :: ncoefs, mass_storage
    REAL(kind=dp), INTENT(OUT), OPTIONAL     :: global_mass
    INTEGER, INTENT(out), OPTIONAL           :: id_nr, ref_count
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_dyn_coeff_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, ncoef_of_kind, nkind
    LOGICAL                                  :: failure
    TYPE(dyn_coeff_type), POINTER            :: coeffs_att

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(dyn_coeff_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(dyn_coeff_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(id_nr)) id_nr = dyn_coeff_set%id_nr
     IF (PRESENT(ref_count)) ref_count = dyn_coeff_set%ref_count
     IF (PRESENT(mass_storage)) mass_storage = dyn_coeff_set%mass_storage
     IF (PRESENT(global_mass)) global_mass = dyn_coeff_set%global_mass
     IF (PRESENT(distribution)) distribution => dyn_coeff_set%distribution
     IF (PRESENT(coeffs_of_kind)) coeffs_of_kind => dyn_coeff_set%coeffs_of_kind
     IF (PRESENT(ncoefs))THEN
       ncoefs=0
       nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
       DO ikind=1, nkind                         
         coeffs_att=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
         IF(ASSOCIATED(coeffs_att))THEN
           CALL get_dyn_coeff(coeffs=coeffs_att,ncoef_of_kind=ncoef_of_kind)
           ncoefs=ncoefs+ncoef_of_kind
         END IF 
       END DO
     END IF
  END IF
END SUBROUTINE get_dyn_coeff_set         

!***************************************************************************
SUBROUTINE dyn_coeff_create(coeffs, error)

    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(coeffs, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_coeff_id=last_coeff_id+1
     coeffs%id_nr=last_coeff_id
     coeffs%ref_count=1
     NULLIFY(coeffs%pos, coeffs%masses, coeffs%vel, coeffs%forces)
  END IF

END SUBROUTINE dyn_coeff_create

!***************************************************************************
SUBROUTINE dyn_coeff_release(coeffs, error)

    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(coeffs)) THEN
     CPPreconditionNoFail(coeffs%ref_count>0,cp_failure_level,routineP,error)
     coeffs%ref_count=coeffs%ref_count-1
     IF (coeffs%ref_count==0) THEN
       IF (ASSOCIATED(coeffs%pos)) THEN 
         DEALLOCATE(coeffs%pos,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(coeffs%masses)) THEN 
         DEALLOCATE(coeffs%masses,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(coeffs%vel)) THEN 
         DEALLOCATE(coeffs%vel,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(coeffs%forces)) THEN 
         DEALLOCATE(coeffs%forces,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       DEALLOCATE(coeffs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(coeffs)

END SUBROUTINE dyn_coeff_release

!***************************************************************************
SUBROUTINE setup_dyn_coeff(coeffs,aux_basis_set,&
                                      n_els, error)

    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    INTEGER, INTENT(in)                      :: n_els
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_dyn_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncgf, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  CALL get_gto_basis_set(gto_basis_set=aux_basis_set,ncgf=ncgf)
  coeffs%ncoef_atom=ncgf
  coeffs%n_els=n_els
  IF(.NOT.ASSOCIATED(coeffs%pos)) THEN          
    ALLOCATE(coeffs%pos(n_els,ncgf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    coeffs%pos(:,:)=1.0_dp
  END IF
  IF(.NOT.ASSOCIATED(coeffs%vel)) THEN          
    ALLOCATE(coeffs%vel(n_els,ncgf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    coeffs%vel(:,:)=0.0_dp
  END IF
  IF(.NOT.ASSOCIATED(coeffs%forces)) THEN          
    ALLOCATE(coeffs%forces(n_els,ncgf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    coeffs%forces(:,:)=0.0_dp
  END IF
   
END SUBROUTINE setup_dyn_coeff
!***************************************************************************
SUBROUTINE get_dyn_coeff(coeffs,pos,vel,forces,masses,n_els,&
                         ncoef_atom,ncoef_of_kind,id_nr,ref_count,error)
   
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pos, vel, forces, masses
    INTEGER, INTENT(OUT), OPTIONAL           :: n_els, ncoef_atom, &
                                                ncoef_of_kind, id_nr, &
                                                ref_count
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_dyn_coeff', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(coeffs),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(id_nr)) id_nr = coeffs%id_nr
     IF (PRESENT(ref_count)) ref_count = coeffs%ref_count
     IF (PRESENT(n_els)) n_els = coeffs%n_els
     IF (PRESENT(ncoef_atom)) ncoef_atom = coeffs%ncoef_atom
     IF (PRESENT(ncoef_of_kind)) ncoef_of_kind = coeffs%ncoef_atom*coeffs%n_els
     IF (PRESENT(pos)) pos => coeffs%pos
     IF (PRESENT(vel)) vel => coeffs%vel
     IF (PRESENT(forces)) forces => coeffs%forces
     IF (PRESENT(masses)) masses => coeffs%masses
  END IF

END SUBROUTINE get_dyn_coeff
!***************************************************************************
SUBROUTINE dyn_coeff_set_initialize(dyn_coeff_set,globenv)         

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER :: group, ikind, n_els, ncoef_atom, ncoef_atom_max, &
      ncoef_of_kind, ncoef_tot, nel_tot, nkind, source, stat,&
      restart_coeff
    REAL(kind=dp)                            :: kg_coeff_target_temp
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: mbuffer
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    type(section_vals_type), pointer         :: coeff_section

!***************************************************************************

  group=globenv%group
  source=globenv%source

! initialize to default values

  dyn_coeff_set%mass_storage=dyn_coeff_mass_global
  dyn_coeff_set%global_mass=100._dp              
  
  nel_tot=0
  ncoef_tot=0
  ncoef_atom_max=0

! get the kg coefficients' parameters from input structure

  NULLIFY(coeff_section)
  coeff_section => section_vals_get_subs_vals(globenv%input_file,"MOTION%MD_NEW%KG_COEFF")
  call section_vals_val_get(coeff_section,"storage",i_val=dyn_coeff_set%mass_storage)
  call section_vals_val_get(coeff_section,"mass",r_val=dyn_coeff_set%global_mass)
  call section_vals_val_get(coeff_section,"restart",i_val=restart_coeff)
  call section_vals_val_get(coeff_section,"temperature",r_val=kg_coeff_target_temp)
    
! get total numbers
  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
  DO ikind=1,nkind
    coeffs=> dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(.not.ASSOCIATED(coeffs)) CYCLE
    CALL get_dyn_coeff(coeffs,n_els=n_els,ncoef_atom=ncoef_atom,&
                       ncoef_of_kind=ncoef_of_kind)
    IF(dyn_coeff_set%mass_storage/=dyn_coeff_mass_global) THEN
      IF(.NOT.ASSOCIATED(coeffs%masses)) THEN
        ALLOCATE(coeffs%masses(n_els,ncoef_atom),STAT=stat)
        IF(stat/=0) CALL stop_memory('routineP','masses',n_els*ncoef_atom)
        coeffs%masses(:,:)=0.0_dp
      END IF
    END IF
    nel_tot=nel_tot+n_els
    ncoef_tot=ncoef_tot+ncoef_of_kind
    ncoef_atom_max=MAX(ncoef_atom_max,ncoef_atom)
  END DO
  CALL mp_sum(nel_tot,group)
  CALL mp_sum(ncoef_tot,group)
  CALL mp_max(ncoef_atom_max,group)

  IF(.not.ALLOCATED(mbuffer)) THEN
    ALLOCATE(mbuffer(nel_tot,ncoef_atom_max),stat=stat) 
    IF(stat/=0) CALL stop_memory('routineP','mbuffer',nel_tot*ncoef_atom_max)
  END IF
  mbuffer(:,:)=0.0_dp 

  dyn_coeff_set%nel_tot=nel_tot 
  dyn_coeff_set%ncoef_tot=ncoef_tot 
  dyn_coeff_set%ncoef_atom_max=ncoef_atom_max

  CALL dyn_coeff_set_init_masses(dyn_coeff_set,mbuffer,globenv)

  CALL dyn_coeff_set_init_vel(kg_coeff_target_temp,dyn_coeff_set,mbuffer,globenv)   

  IF (restart_coeff/=kg_coeff_init_start) THEN
    CALL read_coeff_restart(dyn_coeff_set,globenv,restart_coeff)
  END IF

  DEALLOCATE(mbuffer,stat=stat) 
  IF(stat/=0) CALL stop_memory('routineP','mbuffer')
  
 END SUBROUTINE dyn_coeff_set_initialize
!***************************************************************************
SUBROUTINE read_coeff_restart(dyn_coeff_set,globenv,rtype)          

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    integer, optional                        :: rtype

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coeff_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=200)                       :: file_name
    INTEGER                                  :: group, i, iatom, icoef, iel, &
                                                ikind, istat, n_el, ncoef, &
                                                nkind, nsgf_max, &
                                                restart_unit, source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: cbuffer
    TYPE(dyn_coeff_type), POINTER            :: coeffs

!---------------------------------------------------------------------------

  NULLIFY (coeffs)
  group = globenv%group
  source = globenv%source
  IF (globenv%ionode) THEN
    file_name=globenv%restart_file_name
    CALL open_file(file_name=file_name,&
                   file_action="READ",&
                   file_form="FORMATTED",&
                   file_status="OLD",&
                   unit_number=restart_unit)
    READ (restart_unit,*) n_el, nsgf_max
  END IF
  CALL mp_bcast(n_el,source,group)
  CALL mp_bcast(nsgf_max,source,group)
  ALLOCATE (cbuffer(n_el,nsgf_max),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"cbuffer")
  ALLOCATE (index(n_el),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"index")

  IF (globenv%ionode) THEN
    DO i=1,n_el
      READ (restart_unit,*) index(i),ncoef,cbuffer(i,1:ncoef) 
    END DO
  END IF
  CALL mp_bcast(cbuffer,source,group)
  CALL mp_bcast(index,source,group)
  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
  DO ikind=1, nkind
    coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(ASSOCIATED(coeffs))THEN
      DO iel= 1,coeffs%n_els
        iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
        i=locate(index,iatom)
        DO icoef=1,coeffs%ncoef_atom
          coeffs%pos(iel,icoef)=cbuffer(i,icoef)   
        END DO
      END DO
    END IF
  END DO
  IF((PRESENT(rtype)).AND.(rtype == kg_coeff_pos_vel_restart)) THEN
    IF (globenv%ionode) THEN
      DO i=1,n_el
        READ (restart_unit,*) iatom,ncoef,cbuffer(i,1:ncoef)
      END DO
    END IF
    CALL mp_bcast(cbuffer,source,group)
    DO ikind=1, nkind
      coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF(ASSOCIATED(coeffs))THEN
        DO iel= 1,coeffs%n_els
          iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
          i=locate(index,iatom)
          DO icoef=1,coeffs%ncoef_atom
            coeffs%vel(iel,icoef)=cbuffer(i,icoef)  
          END DO
        END DO
      END IF
    END DO
  END IF
  DEALLOCATE(cbuffer,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"cbuffer")
  DEALLOCATE(index,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"index")

  IF (globenv%ionode) THEN
   CALL close_file(unit_number=restart_unit)
  END IF

  END SUBROUTINE read_coeff_restart
!***************************************************************************
SUBROUTINE write_coeff_restart(dyn_coeff_set,atomic_kind_set,globenv)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'write_coeff_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=200)                       :: file_name
    INTEGER :: group, i, iatom, icoef, iel, ikind, ioff, istat, n_el, natom, &
      nc, nkind, nsgf, nsgf_max, offset, restart_unit, source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, ncoef
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: cbuffer
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set

!---------------------------------------------------------------------------

  NULLIFY (coeffs,aux_basis_set,atomic_kind)
  group = globenv%group
  source = globenv%source
  nsgf_max=0
  n_el=0
  IF (globenv%ionode) THEN
    file_name=globenv%restart_file_name
    CALL open_file(file_name=file_name,&
                   file_action="WRITE",&
                   file_form="FORMATTED",&
                   file_status="REPLACE",&
                   unit_number=restart_unit)
    DO ikind=1,size(atomic_kind_set)
      atomic_kind=> atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           aux_basis_set=aux_basis_set,&
                           natom=natom)
      IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                             nsgf=nsgf)
      nsgf_max=MAX(nsgf_max,nsgf)
      n_el = n_el + natom
    END DO
    WRITE (restart_unit,*) n_el, nsgf_max
  END IF

  CALL mp_bcast(n_el,source,group)
  CALL mp_bcast(nsgf_max,source,group)
  ALLOCATE (cbuffer(n_el,nsgf_max),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"cbuffer")
  ALLOCATE (index(n_el),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"index")
  ALLOCATE (ncoef(n_el),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"ncoef")
  cbuffer(:,:)=0._dp
  index(:)=0
  ncoef(:)=0

  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
  offset = 0
  DO ikind=1, nkind
    coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    atomic_kind=> atomic_kind_set(ikind)
    IF(.NOT.ASSOCIATED(coeffs))THEN
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom)
       offset = offset + natom
    ELSE
      DO iel= 1,coeffs%n_els
        iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
        ioff=iatom - offset
        index(ioff)=iatom
        ncoef(ioff)=coeffs%ncoef_atom
        DO icoef=1,coeffs%ncoef_atom
          cbuffer(ioff,icoef)=coeffs%pos(iel,icoef)   
        END DO
      END DO
    END IF
  END DO
  CALL mp_sum(cbuffer,group)
  CALL mp_sum(index,group)
  CALL mp_sum(ncoef,group)
  IF (globenv%ionode) THEN
    DO i=1,n_el
      nc=ncoef(i)
      WRITE(restart_unit,'(2I6,6f20.10)') index(i),ncoef(i),cbuffer(i,1:nc)
    END DO
  END IF
  DEALLOCATE(cbuffer,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"cbuffer")
  DEALLOCATE(index,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"index")
  DEALLOCATE(ncoef,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"ncoef")

  IF (globenv%ionode) THEN
   CALL close_file(unit_number=restart_unit)
  END IF

  END SUBROUTINE write_coeff_restart

!***************************************************************************
SUBROUTINE dyn_coeff_set_dump(dyn_coeff_set,atomic_kind_set,error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: file_name = 'coeffs_dump', &
      routineN = 'dyn_coeff_set_dump', routineP = moduleN//':'//routineN

    INTEGER                                  :: cdump_unit, group, i, iatom, &
                                                icoef, iel, ikind, ioff, &
                                                istat, natom, nc, nkind, &
                                                offset
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, ncoef
    LOGICAL                                  :: should_write
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: cbuffer, fbuffer, vbuffer
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dyn_coeff_type), POINTER            :: coeffs

!---------------------------------------------------------------------------

  NULLIFY (coeffs,atomic_kind)
  logger => cp_error_get_logger(error)
  should_write=(logger%para_env%mepos==logger%para_env%source)
  group = logger%para_env%group
  IF (should_write) THEN
    CALL open_file(file_name=file_name,&
                   file_action="WRITE",&
                   file_form="FORMATTED",&
                   file_status="REPLACE",&
                   unit_number=cdump_unit)
    WRITE (cdump_unit,*) dyn_coeff_set%nel_tot, dyn_coeff_set%ncoef_atom_max
  END IF

  ALLOCATE (cbuffer(dyn_coeff_set%nel_tot,dyn_coeff_set%ncoef_atom_max),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"cbuffer")
  ALLOCATE (vbuffer(dyn_coeff_set%nel_tot,dyn_coeff_set%ncoef_atom_max),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"vbuffer")
  ALLOCATE (fbuffer(dyn_coeff_set%nel_tot,dyn_coeff_set%ncoef_atom_max),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"fbuffer")
  ALLOCATE (index(dyn_coeff_set%nel_tot),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"index")
  ALLOCATE (ncoef(dyn_coeff_set%nel_tot),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"ncoef")
  cbuffer(:,:)=0._dp
  vbuffer(:,:)=0._dp
  fbuffer(:,:)=0._dp
  index(:)=0
  ncoef(:)=0

  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
  offset = 0
  DO ikind=1, nkind
    coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    atomic_kind=> atomic_kind_set(ikind)
    IF(.NOT.ASSOCIATED(coeffs))THEN
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom)
       offset = offset + natom
    ELSE
      DO iel= 1,coeffs%n_els
        iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
        ioff=iatom - offset
        index(ioff)=iatom
        ncoef(ioff)=coeffs%ncoef_atom
        DO icoef=1,coeffs%ncoef_atom
          cbuffer(ioff,icoef)=coeffs%pos(iel,icoef)   
          vbuffer(ioff,icoef)=coeffs%vel(iel,icoef)   
          fbuffer(ioff,icoef)=coeffs%forces(iel,icoef)   
        END DO
      END DO
    END IF
  END DO
  CALL mp_sum(cbuffer,group)
  CALL mp_sum(vbuffer,group)
  CALL mp_sum(fbuffer,group)
  CALL mp_sum(index,group)
  CALL mp_sum(ncoef,group)
  IF (should_write) THEN
    DO i=1,dyn_coeff_set%nel_tot
      nc=dyn_coeff_set%ncoef_atom_max
      WRITE(cdump_unit,'(2I6,6f20.10)') index(i),ncoef(i),cbuffer(i,1:nc)
    END DO
    DO i=1,dyn_coeff_set%nel_tot
      nc=dyn_coeff_set%ncoef_atom_max
      WRITE(cdump_unit,'(2I6,6f20.10)') index(i),ncoef(i),vbuffer(i,1:nc)
    END DO
    DO i=1,dyn_coeff_set%nel_tot
      nc=dyn_coeff_set%ncoef_atom_max
      WRITE(cdump_unit,'(2I6,6f20.10)') index(i),ncoef(i),fbuffer(i,1:nc)
    END DO
  END IF
  DEALLOCATE(cbuffer,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"cbuffer")
  DEALLOCATE(vbuffer,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"vbuffer")
  DEALLOCATE(fbuffer,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"fbuffer")
  DEALLOCATE(index,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"index")
  DEALLOCATE(ncoef,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"ncoef")

  IF (should_write) THEN
   CALL close_file(unit_number=cdump_unit)
  END IF

  END SUBROUTINE dyn_coeff_set_dump

!!*************************************************************************
SUBROUTINE dyn_coeff_set_init_masses(dyn_coeff_set,mbuffer,globenv)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: mbuffer
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: file_name = 'CMASSES', &
      routineN = 'dyn_coeff_set_init_masses', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cmass_unit, group, i, iatom, &
                                                icoef, iel, ifin, ikind, &
                                                ioff, istart, nkind, source, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: ckind, nel_kind
    TYPE(dyn_coeff_type), POINTER            :: coeffs

  group=globenv%group
  source=globenv%source
  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)

  ALLOCATE (index(dyn_coeff_set%nel_tot),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineP,"index")
  index(:)=0

! initialize coefficient effective masses


  SELECT CASE (dyn_coeff_set%mass_storage)
  CASE (dyn_coeff_mass_kind)
    IF(.not.ALLOCATED(ckind)) THEN
      ALLOCATE(ckind(nkind),stat=stat) 
      IF(stat/=0) CALL stop_memory('routineP','ckind',nkind)
    END IF
    IF(.not.ALLOCATED(nel_kind)) THEN
      ALLOCATE(nel_kind(nkind),stat=stat) 
      IF(stat/=0) CALL stop_memory('routineP','nel_kind',nkind)
    END IF
    IF (globenv%ionode) THEN
      CALL open_file(file_name=file_name,&
                   file_action="READ",&
                   file_form="FORMATTED",&
                   unit_number=cmass_unit)
      READ (cmass_unit,*) ckind(1:nkind)
      READ (cmass_unit,*) nel_kind(1:nkind)
    END IF
    CALL mp_bcast(ckind,source,group)
    DO ikind=1,nkind
      coeffs=> dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF(.not.ASSOCIATED(coeffs)) CYCLE
      coeffs%masses(:,:)=ckind(ikind)
    END DO
    IF (globenv%ionode) THEN
      ioff=0
      DO ikind=1,nkind
        istart=ioff+1
        ifin=ioff+nel_kind(ikind)
        coeffs=> dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
        IF(ASSOCIATED(coeffs)) THEN
          mbuffer(istart:ifin,:)=ckind(ikind)  
          ioff =ioff+nel_kind(ikind)
        END IF
      END DO
    END IF
    DEALLOCATE(ckind,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"ckind")
    DEALLOCATE(nel_kind,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"nel_kind")
  CASE (dyn_coeff_mass_atom)
    IF (globenv%ionode) THEN
      CALL open_file(file_name=file_name,&
                   file_action="READ",&
                   file_form="FORMATTED",&
                   unit_number=cmass_unit)
      DO i=1,dyn_coeff_set%nel_tot 
        READ (cmass_unit,*) index(i),mbuffer(i,1:1)
      END DO
    END IF
    CALL mp_bcast(mbuffer,source,group)
    CALL mp_bcast(index,source,group)
    DO ikind=1, nkind
      coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF(ASSOCIATED(coeffs))THEN
        DO iel= 1,coeffs%n_els
          iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
          i=locate(index,iatom)
          DO icoef=1,coeffs%ncoef_atom
            coeffs%masses(iel,icoef)=mbuffer(i,1)
            mbuffer(i,icoef)=mbuffer(i,1)
          END DO
        END DO
      END IF
    END DO
    IF (globenv%ionode) THEN
      CALL close_file(unit_number=cmass_unit)
    END IF
  CASE (dyn_coeff_mass_coeff)
    IF (globenv%ionode) THEN
      CALL open_file(file_name=file_name,&
                   file_action="READ",&
                   file_form="FORMATTED",&
                   unit_number=cmass_unit)
      DO i=1,dyn_coeff_set%nel_tot
        READ (cmass_unit,*) index(i),mbuffer(i,1:dyn_coeff_set%ncoef_atom_max)
      END DO
    END IF
    CALL mp_bcast(mbuffer,source,group)
    CALL mp_bcast(index,source,group)
    DO ikind=1, nkind
      coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF(ASSOCIATED(coeffs))THEN
        DO iel= 1,coeffs%n_els
          iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
          i=locate(index,iatom)
          DO icoef=1,coeffs%ncoef_atom
             coeffs%masses(iel,icoef)=mbuffer(i,icoef)
          END DO
        END DO
      END IF
    END DO
    IF (globenv%ionode) THEN
      CALL close_file(unit_number=cmass_unit)
    END IF
  END SELECT

  DEALLOCATE (index,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineP,"index")

END SUBROUTINE dyn_coeff_set_init_masses
!*************************************************************************
SUBROUTINE dyn_coeff_set_init_vel(kg_coeff_target_temp,dyn_coeff_set,mbuffer,globenv)

    REAL(KIND=dp)                            :: kg_coeff_target_temp
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: mbuffer
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_init_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, iatom, icoef, iel, &
                                                ikind, ioff, istat, nkind, &
                                                offset, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: n_els_of_kind
    LOGICAL                                  :: gmass
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: vbuffer
    REAL(KIND=dp)                            :: denom, ekin, scale, &
                                                temp_coefs, vcom 
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(unit_convert_type), POINTER         :: unit

  NULLIFY (coeffs)
  gmass=(dyn_coeff_set%mass_storage==dyn_coeff_mass_global)
  group = globenv%group

  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)      
  vcom=0._dp
  denom=0._dp

  ALLOCATE (vbuffer(dyn_coeff_set%nel_tot,dyn_coeff_set%ncoef_atom_max),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"vbuffer")
  vbuffer(:,:)=0._dp

  ALLOCATE (n_els_of_kind(nkind),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"n_els_of_kind")
  n_els_of_kind(:)=0

! calculate com velocities

  IF (globenv%ionode) THEN
    DO iel=1,dyn_coeff_set%nel_tot
      DO icoef=1,dyn_coeff_set%ncoef_atom_max
        vbuffer(iel,icoef)=gasdev(globenv%idum)
        IF (gmass) THEN
          vcom=vcom+dyn_coeff_set%global_mass*vbuffer(iel,icoef)
          denom=denom+dyn_coeff_set%global_mass
        ELSE
          vcom=vcom+mbuffer(iel,icoef)*vbuffer(iel,icoef)
          denom=denom+mbuffer(iel,icoef)
        END IF
      END DO
    END DO
    vcom=vcom/denom
    ekin=0._dp

    DO iel=1,dyn_coeff_set%nel_tot 
      DO icoef=1,dyn_coeff_set%ncoef_atom_max
        vbuffer(iel,icoef)=vbuffer(iel,icoef)-vcom
        IF (gmass) THEN
          ekin=ekin+0.5_dp*dyn_coeff_set%global_mass*vbuffer(iel,icoef) &
               *vbuffer(iel,icoef)
        ELSE
          ekin=ekin+0.5_dp*mbuffer(iel,icoef)*vbuffer(iel,icoef)*vbuffer(iel,icoef)
        END IF
      END DO
    END DO

! calculate coeff temperature
    temp_coefs = 2.0_dp * ekin / REAL ( dyn_coeff_set%ncoef_tot, dp )

! scale velocties to get the correct initial coefficent temperature
 
    scale= SQRT ( kg_coeff_target_temp / temp_coefs )
    DO iel=1,dyn_coeff_set%nel_tot
      DO icoef=1,dyn_coeff_set%ncoef_atom_max
           vbuffer(iel,icoef)=vbuffer(iel,icoef)*scale
      END DO
    END DO

! get the initial coef temperature
    ekin=0._dp
    vcom=0._dp
    DO iel=1,dyn_coeff_set%nel_tot
      DO icoef=1,dyn_coeff_set%ncoef_atom_max
        IF (gmass) THEN
          vcom=vcom+dyn_coeff_set%global_mass*vbuffer(iel,icoef)
          ekin=ekin+0.5_dp*dyn_coeff_set%global_mass*vbuffer(iel,icoef) &
               *vbuffer(iel,icoef)
        ELSE
          vcom=vcom+mbuffer(iel,icoef)*vbuffer(iel,icoef)
          ekin=ekin+0.5_dp*mbuffer(iel,icoef)*vbuffer(iel,icoef)*vbuffer(iel,icoef)
        END IF
      END DO
    END DO

  END IF
  CALL mp_bcast(vbuffer,globenv%source,group)
  vcom=vcom/denom

  DO ikind=1, nkind
    n_els_of_kind(ikind) = dyn_coeff_set%distribution%n_el(ikind)
  END DO
  CALL mp_sum(n_els_of_kind,group)
   
  offset=0
  DO ikind=1, nkind
    coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(.NOT.ASSOCIATED(coeffs))THEN
       offset = offset + n_els_of_kind(ikind)
    ELSE
      DO iel= 1,coeffs%n_els
        iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
        ioff=iatom - offset
        DO icoef=1,coeffs%ncoef_atom
          coeffs%vel(iel,icoef)=vbuffer(ioff,icoef)
        END DO
      END DO
    END IF
  END DO

  DEALLOCATE(vbuffer,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"vbuffer")
  DEALLOCATE(n_els_of_kind,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineP,"n_els_of_kind")
   
  IF ( globenv % ionode .AND. globenv % print_level > 1 ) THEN
     CALL get_cp2k_units ( unit )
     temp_coefs = 2.0_dp * ekin / REAL ( dyn_coeff_set%ncoef_tot, dp )
     WRITE ( globenv%scr, '( A, T61, F18.2, A2 )' ) &
          ' Initial Coef Temperature ', temp_coefs / unit % kelvin, " K"
     WRITE ( globenv%scr, '( A, T61, F20.12 )' ) &
          'Coefs Centre of mass velocity ', vcom
  END IF

END SUBROUTINE dyn_coeff_set_init_vel
!*****************************************************************************

END MODULE dynamical_coeff_types
