!---------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
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
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dyn_coeff_mass_atom,&
                                             dyn_coeff_mass_coeff,&
                                             dyn_coeff_mass_global,&
                                             dyn_coeff_mass_kind
  USE input_section_types,             ONLY: section_get_cval,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: gasdev,&
                                             locate
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dynamical_coeff_types'
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
 
  PUBLIC :: dyn_coeff_set_create, dyn_coeff_set_initialize,&
            dyn_coeff_set_release, dyn_coeff_set_retain,&
            get_dyn_coeff_set, dyn_coeff_create, setup_dyn_coeff, &
            dyn_coeff_release, get_dyn_coeff


! *****************************************************************************
  CONTAINS

! *****************************************************************************

!!*** **********************************************************************
SUBROUTINE dyn_coeff_set_create(dyn_coeff_set,atomic_kind_set,&
                                distribution, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
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
     nkind=SIZE(atomic_kind_set)
     IF (PRESENT(distribution)) THEN
       IF (ASSOCIATED(distribution)) THEN
         ALLOCATE(dyn_coeff_set%coeffs_of_kind(nkind), stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         dyn_coeff_set%distribution_method=dyn_coeff_distributed
         DO ikind=1,nkind 
            NULLIFY(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind,aux_basis_set=aux_basis_set)
            n_els = distribution%n_el(ikind)
            IF ((.NOT.ASSOCIATED(aux_basis_set)).OR.(n_els ==0)) CYCLE
!dbg           write(6,*) 'dyncoeffsetcreate: n_els, ikind', n_els, ikind 
            CALL dyn_coeff_create(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                             error=error)
            CALL setup_dyn_coeff(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs,&
               aux_basis_set=aux_basis_set, n_els=n_els, error=error)
         END DO
         CALL distribution_1d_retain(distribution,error=error)
         CALL distribution_1d_release(dyn_coeff_set%distribution,error=error)
         dyn_coeff_set%distribution => distribution
       END IF
     ELSE 
       dyn_coeff_set%distribution_method=dyn_coeff_replicated
       CALL stop_program ('dyn_coeff_set_create','replicated coefs NYI')
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
SUBROUTINE dyn_coeff_set_initialize(dyn_coeff_set,globenv,work_section,skipvel,error)         

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: work_section
    LOGICAL, INTENT(IN), OPTIONAL            :: skipvel
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER :: group, i, iatom, icoef, iel, ikind, n_els, ncoef_atom, &
      ncoef_atom_max, ncoef_of_kind, ncoef_tot, nel_tot, nkind, source, stat
    INTEGER, DIMENSION(:), POINTER           :: ibuffer
    LOGICAL                                  :: explicit, my_skipvel, restart
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: mbuffer
    REAL(kind=dp)                            :: kg_coeff_target_temp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(section_vals_type), POINTER         :: coeff_section, &
                                                storage_section, work_section2

!***************************************************************************

  group=globenv%group
  source=globenv%source
  restart = .FALSE.
  my_skipvel = .FALSE.
  IF (PRESENT(skipvel)) my_skipvel = skipvel
  NULLIFY(work_section2, buffer,ibuffer)
  ! initialize to default values

  dyn_coeff_set%mass_storage=dyn_coeff_mass_global
  dyn_coeff_set%global_mass=100._dp              

  nel_tot=0
  ncoef_tot=0
  ncoef_atom_max=0

  ! get the kg coefficients' parameters from input structure

  NULLIFY(coeff_section)
  coeff_section => section_vals_get_subs_vals(globenv%input_file,"MOTION%MD%KG_COEFF")
  storage_section =>  section_vals_get_subs_vals(coeff_section,"STORAGE")
  CALL section_vals_val_get(storage_section,"type",i_val=dyn_coeff_set%mass_storage)
  CALL section_vals_val_get(coeff_section,"mass",r_val=dyn_coeff_set%global_mass)
  CALL section_vals_val_get(coeff_section,"temperature",r_val=kg_coeff_target_temp)

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

  dyn_coeff_set%nel_tot=nel_tot 
  dyn_coeff_set%ncoef_tot=ncoef_tot 
  dyn_coeff_set%ncoef_atom_max=ncoef_atom_max

  ALLOCATE(mbuffer(nel_tot,ncoef_atom_max),stat=stat) 
  IF(stat/=0) CALL stop_memory('routineP','mbuffer',nel_tot*ncoef_atom_max)
  mbuffer(:,:)=0.0_dp 
  CALL dyn_coeff_set_init_masses(dyn_coeff_set,mbuffer,storage_section,globenv,error)


  work_section2 => section_vals_get_subs_vals(work_section,"COORD",error=error)
  CALL section_vals_get(work_section2, explicit=explicit, error=error)
  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
  IF (explicit) THEN
     CALL section_vals_val_get(work_section,"RESTART_INDEX%_DEFAULT_KEYWORD_",i_vals=ibuffer,error=error)
     DO ikind=1, nkind
        coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
        IF(ASSOCIATED(coeffs))THEN
           DO iel= 1,coeffs%n_els
              iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
              i=locate(ibuffer,iatom)
              CALL section_vals_val_get(work_section,"COORD%_DEFAULT_KEYWORD_",r_vals=buffer,&
                   i_rep_val=i,error=error)
              DO icoef=1,coeffs%ncoef_atom
                 coeffs%pos(iel,icoef)=buffer(icoef)   
              END DO
           END DO
        END IF
     END DO
  END IF

  IF (.NOT.my_skipvel) THEN
     work_section2 => section_vals_get_subs_vals(work_section,"VELOCITY",error=error)
     CALL section_vals_get(work_section2, explicit=explicit, error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(work_section,"RESTART_INDEX%_DEFAULT_KEYWORD_",i_vals=ibuffer,error=error)
        DO ikind=1, nkind
           coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
           IF(ASSOCIATED(coeffs))THEN
              DO iel= 1,coeffs%n_els
                 iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
                 i=locate(ibuffer,iatom)
                 CALL section_vals_val_get(work_section,"VELOCITY%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
                 DO icoef=1,coeffs%ncoef_atom
                    coeffs%vel(iel,icoef)=buffer(icoef)  
                 END DO
              END DO
           END IF
        END DO
     ELSE
        CALL dyn_coeff_set_init_vel(kg_coeff_target_temp,dyn_coeff_set,mbuffer,globenv)   
     END IF
  END IF

  DEALLOCATE(mbuffer,stat=stat) 
  IF(stat/=0) CALL stop_memory('routineP','mbuffer')

END SUBROUTINE dyn_coeff_set_initialize


!!*************************************************************************
  SUBROUTINE dyn_coeff_set_init_masses(dyn_coeff_set,mbuffer,storage_section,globenv,error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: mbuffer
    TYPE(section_vals_type), POINTER         :: storage_section
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dyn_coeff_set_init_masses', &
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

    IF (dyn_coeff_set%mass_storage /= dyn_coeff_mass_global) THEN
       ALLOCATE (INDEX(dyn_coeff_set%nel_tot),STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineP,"index")
       INDEX(:)=0

       IF (globenv%ionode) THEN
          CALL open_file(file_name=section_get_cval(storage_section,"FILENAME",error=error),&
               file_action="READ",&
               file_form="FORMATTED",&
               unit_number=cmass_unit)
       END IF
       ! initialize coefficient effective masses
       SELECT CASE (dyn_coeff_set%mass_storage)
       CASE (dyn_coeff_mass_kind)
          IF(.NOT.ALLOCATED(ckind)) THEN
             ALLOCATE(ckind(nkind),stat=stat) 
             IF(stat/=0) CALL stop_memory('routineP','ckind',nkind)
          END IF
          IF(.NOT.ALLOCATED(nel_kind)) THEN
             ALLOCATE(nel_kind(nkind),stat=stat) 
             IF(stat/=0) CALL stop_memory('routineP','nel_kind',nkind)
          END IF
          IF (globenv%ionode) THEN
             READ (cmass_unit,*) ckind(1:nkind)
             READ (cmass_unit,*) nel_kind(1:nkind)
          END IF
          CALL mp_bcast(ckind,source,group)
          DO ikind=1,nkind
             coeffs=> dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
             IF(.NOT.ASSOCIATED(coeffs)) CYCLE
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
             DO i=1,dyn_coeff_set%nel_tot 
                READ (cmass_unit,*) INDEX(i),mbuffer(i,1:1)
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
       CASE (dyn_coeff_mass_coeff)
          IF (globenv%ionode) THEN
             DO i=1,dyn_coeff_set%nel_tot
                READ (cmass_unit,*) INDEX(i),mbuffer(i,1:dyn_coeff_set%ncoef_atom_max)
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
       END SELECT
       IF (globenv%ionode) THEN
          CALL close_file(unit_number=cmass_unit)
       END IF
       DEALLOCATE (index,STAT=stat)
       IF (stat /= 0) CALL stop_memory(routineP,"index")
    END IF
  END SUBROUTINE dyn_coeff_set_init_masses
!*************************************************************************
SUBROUTINE dyn_coeff_set_init_vel(kg_coeff_target_temp,dyn_coeff_set,mbuffer,globenv)

    REAL(KIND=dp)                            :: kg_coeff_target_temp
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: mbuffer
    TYPE(global_environment_type), POINTER   :: globenv

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

  NULLIFY (coeffs,unit)
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
    vcom=vcom/denom

  END IF
  CALL mp_bcast(vbuffer,globenv%source,group)

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
