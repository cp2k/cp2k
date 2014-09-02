!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular auxamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the data structure for the auxilliary basis SCP coefficents .
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
MODULE scp_coeff_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE f77_blas
  USE input_constants,                 ONLY: copy_coeff,&
                                             copy_fcoeff,&
                                             use_aux_basis_set
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_coeff_types'
  INTEGER, PRIVATE, SAVE :: last_coeff_id=0
  INTEGER, PRIVATE, SAVE :: last_coeff_set_id=0

  INTEGER, PARAMETER,PUBLIC :: aux_coeff_replicated=1, aux_coeff_distributed=2

! *****************************************************************************
  TYPE aux_coeff_set_type
   INTEGER :: ref_count, id_nr
   INTEGER :: distribution_method
   INTEGER :: ncoef_atom_max, ncoef_tot, nel_tot
   LOGICAL :: propagate
   TYPE ( distribution_1d_type ), POINTER :: distribution
   TYPE ( aux_coeff_p_type ), POINTER, DIMENSION(:) :: coeffs_of_kind
  END TYPE aux_coeff_set_type

! *****************************************************************************
  TYPE aux_coeff_set_p_type
   TYPE(aux_coeff_set_type), POINTER :: aux_coeff_set
  END TYPE aux_coeff_set_p_type

! *****************************************************************************
  TYPE aux_coeff_type
   INTEGER :: ref_count, id_nr
   INTEGER :: n_els
   INTEGER :: ncoef_atom
   REAL(KIND=dp), DIMENSION(:,:), POINTER :: c, fc
  END TYPE aux_coeff_type

! *****************************************************************************
  TYPE aux_coeff_p_type
   TYPE(aux_coeff_type), POINTER :: coeffs
  END TYPE aux_coeff_p_type

! *** Public structures ***

  PUBLIC :: aux_coeff_set_type,&
            aux_coeff_type

! *** Public Subroutines ***

  PUBLIC :: aux_coeff_set_release, aux_coeff_set_retain,&
            get_aux_coeff

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param aux_coeff_set ...
!> \param qs_kind_set ...
!> \param distribution ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_set_create(aux_coeff_set,qs_kind_set,&
       distribution, error)

    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: distribution
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_set_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, n_els, nkind, stat
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    failure=.FALSE.

    NULLIFY(qs_kind,aux_basis_set)
    ALLOCATE(aux_coeff_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_coeff_set_id=last_coeff_set_id+1
       aux_coeff_set%id_nr=last_coeff_set_id
       aux_coeff_set%ref_count=1
       NULLIFY(aux_coeff_set%distribution, aux_coeff_set%coeffs_of_kind)
       !    set coeffs_of_kind
       nkind=SIZE(qs_kind_set)
       IF (PRESENT(distribution)) THEN
          IF (ASSOCIATED(distribution)) THEN
             ALLOCATE(aux_coeff_set%coeffs_of_kind(nkind), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             aux_coeff_set%distribution_method=aux_coeff_distributed
             DO ikind=1,nkind
                NULLIFY(aux_coeff_set%coeffs_of_kind(ikind)%coeffs)
                qs_kind => qs_kind_set(ikind)
                CALL get_qs_kind(qs_kind,aux_basis_set=aux_basis_set)
                n_els = distribution%n_el(ikind)
                IF ((.NOT.ASSOCIATED(aux_basis_set)).OR.(n_els ==0)) CYCLE
                CALL aux_coeff_create(aux_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                     error=error)
                CALL setup_aux_coeff(aux_coeff_set%coeffs_of_kind(ikind)%coeffs,&
                     aux_basis_set=aux_basis_set, n_els=n_els, error=error)
             END DO

             CALL distribution_1d_retain(distribution,error=error)
             CALL distribution_1d_release(aux_coeff_set%distribution,error=error)
             aux_coeff_set%distribution => distribution
          END IF
       ELSE
          aux_coeff_set%distribution_method=aux_coeff_replicated
          CALL stop_program(routineN,moduleN,__LINE__,'Replicated coeffs NYI')
       END IF
    END IF

  END SUBROUTINE aux_coeff_set_create

! *****************************************************************************
!> \brief ...
!> \param aux_coeff_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_set_retain(aux_coeff_set, error)

    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_set_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(aux_coeff_set),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(aux_coeff_set%ref_count>0,cp_failure_level,routineP,error)
       aux_coeff_set%ref_count=aux_coeff_set%ref_count+1
    END IF

  END SUBROUTINE aux_coeff_set_retain

! *****************************************************************************
!> \brief ...
!> \param aux_coeff_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_zero_fc (aux_coeff_set, error)

    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_zero_fc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: force
    TYPE(aux_coeff_type), POINTER            :: local_coeffs

    failure=.FALSE.
    NULLIFY ( local_coeffs, force )

    CPPrecondition(ASSOCIATED(aux_coeff_set),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO ikind=1,SIZE(aux_coeff_set%coeffs_of_kind)
          local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
          IF ( .NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
          CALL get_aux_coeff(coeffs=local_coeffs,fc=force,error=error)
          force(:,:)=0.0_dp
       END DO
    END IF

  END SUBROUTINE aux_coeff_zero_fc

! *****************************************************************************
!> \brief ...
!> \param aux_coeff_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_set_release(aux_coeff_set, error)

    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_set_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(aux_coeff_set)) THEN
       CPPreconditionNoFail(aux_coeff_set%ref_count>0,cp_failure_level,routineP,error)
       aux_coeff_set%ref_count=aux_coeff_set%ref_count-1
       IF (aux_coeff_set%ref_count==0) THEN
          CALL distribution_1d_release(aux_coeff_set%distribution, error=error)
          IF(ASSOCIATED(aux_coeff_set%coeffs_of_kind)) THEN
             DO i=1,SIZE(aux_coeff_set%coeffs_of_kind)
                IF(ASSOCIATED(aux_coeff_set%coeffs_of_kind(i)%coeffs)) THEN
                   CALL aux_coeff_release(aux_coeff_set%coeffs_of_kind(i)%coeffs,error=error)
                ENDIF
             END DO
             DEALLOCATE(aux_coeff_set%coeffs_of_kind,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(aux_coeff_set, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(aux_coeff_set)

  END SUBROUTINE aux_coeff_set_release

! *****************************************************************************
!> \brief ...
!> \param aux_coeff_set ...
!> \param distribution ...
!> \param coeffs_of_kind ...
!> \param ncoefs ...
!> \param id_nr ...
!> \param ref_count ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_aux_coeff_set(aux_coeff_set,distribution, coeffs_of_kind,&
       ncoefs,id_nr,ref_count,error)

    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: distribution
    TYPE(aux_coeff_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: coeffs_of_kind
    INTEGER, INTENT(out), OPTIONAL           :: ncoefs, id_nr, ref_count
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_aux_coeff_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, ncoef_of_kind, nkind
    LOGICAL                                  :: failure
    TYPE(aux_coeff_type), POINTER            :: coeffs_att

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(aux_coeff_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(aux_coeff_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(id_nr)) id_nr = aux_coeff_set%id_nr
       IF (PRESENT(ref_count)) ref_count = aux_coeff_set%ref_count
       IF (PRESENT(distribution)) distribution => aux_coeff_set%distribution
       IF (PRESENT(coeffs_of_kind)) coeffs_of_kind => aux_coeff_set%coeffs_of_kind
       IF (PRESENT(ncoefs))THEN
          ncoefs=0
          nkind=SIZE(aux_coeff_set%coeffs_of_kind)
          DO ikind=1, nkind
             coeffs_att=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
             IF(ASSOCIATED(coeffs_att))THEN
                CALL get_aux_coeff(coeffs=coeffs_att,ncoef_of_kind=ncoef_of_kind,error=error)
                ncoefs=ncoefs+ncoef_of_kind
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE get_aux_coeff_set

! *****************************************************************************
!> \brief ...
!> \param coeffs ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_create(coeffs, error)

    TYPE(aux_coeff_type), POINTER            :: coeffs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_create', &
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
       NULLIFY( coeffs%c, coeffs%fc )
    END IF

  END SUBROUTINE aux_coeff_create

! *****************************************************************************
!> \brief ...
!> \param coeffs ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_release(coeffs, error)

    TYPE(aux_coeff_type), POINTER            :: coeffs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(coeffs)) THEN
       CPPreconditionNoFail(coeffs%ref_count>0,cp_failure_level,routineP,error)
       coeffs%ref_count=coeffs%ref_count-1
       IF (coeffs%ref_count==0) THEN
          IF (ASSOCIATED(coeffs%c)) THEN
             DEALLOCATE(coeffs%c,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(coeffs%fc)) THEN
             DEALLOCATE(coeffs%fc,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(coeffs,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(coeffs)

  END SUBROUTINE aux_coeff_release

! *****************************************************************************
!> \brief ...
!> \param coeffs ...
!> \param aux_basis_set ...
!> \param n_els ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE setup_aux_coeff(coeffs,aux_basis_set, n_els, error)

    TYPE(aux_coeff_type), POINTER            :: coeffs
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    INTEGER, INTENT(in)                      :: n_els
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_aux_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nsgf, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL get_gto_basis_set(gto_basis_set=aux_basis_set,nsgf=nsgf)
    coeffs%ncoef_atom=nsgf
    coeffs%n_els=n_els
    IF(.NOT.ASSOCIATED(coeffs%c)) THEN
       ALLOCATE(coeffs%c(n_els,nsgf),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       coeffs%c(:,:)=0.0_dp
    END IF
    IF(.NOT.ASSOCIATED(coeffs%fc)) THEN
       ALLOCATE(coeffs%fc(n_els,nsgf),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       coeffs%fc(:,:)=0.0_dp
    END IF

  END SUBROUTINE setup_aux_coeff
! *****************************************************************************
!> \brief ...
!> \param rep_coeff ...
!> \param qs_kind_set ...
!> \param local_particles ...
!> \param aux_coeff_set ...
!> \param para_env ...
!> \param copy_type_id ...
!> \param copy_in ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_replicated_coeff ( rep_coeff, atomic_kind_set, qs_kind_set, &
                                      local_particles, aux_coeff_set, para_env, &
                                      copy_type_id, copy_in, error )

    REAL(dp), DIMENSION(:, :), POINTER       :: rep_coeff
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(qs_kind_type), POINTER              :: qs_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: copy_type_id
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_in
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_replicated_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ii, ikind, &
                                                iparticle_local, maxgtops, &
                                                natoms, nparticle_local, &
                                                stat, ub
    LOGICAL                                  :: failure, my_copy_in
    REAL(dp), POINTER                        :: c( :, : ), fc( :, : )
    TYPE(aux_coeff_type), POINTER            :: local_coeffs

  NULLIFY ( c, fc )
  my_copy_in = .TRUE.
  IF ( PRESENT ( copy_in ) ) my_copy_in = copy_in

! Ensuring a blank canvas
  IF ( ASSOCIATED ( rep_coeff ) ) THEN
   CALL deallocate_replicated_coeff ( rep_coeff, qs_kind_set, local_particles, &
                                      aux_coeff_set, copy_type_id, .FALSE., error )

  END IF

! Getting dimensions of the replicated coeff array
  CALL get_atomic_kind_set(atomic_kind_set, natom=natoms)
  CALL get_qs_kind_set ( qs_kind_set, maxgtops=maxgtops, basis_set_id=use_aux_basis_set)

! Allocating the pointer: MAXGTOPS is the maximum number of spherical orbitals in a set * number of sets
  ALLOCATE(rep_coeff ( natoms, maxgtops ), stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineN,error,failure)
! Zero
  rep_coeff = 0._dp
  IF ( my_copy_in ) THEN
    DO ikind = 1, SIZE ( qs_kind_set )
      local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF (.NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
  ! Get the number of particles of ikind ( local )
      nparticle_local = local_particles%n_el(ikind)
  ! Loop over the (local) particles
      DO iparticle_local=1,nparticle_local
         iatom = local_particles%list(ikind)%array(iparticle_local)
  ! Assigning the local coefficients:
         SELECT CASE ( copy_type_id )
         CASE ( copy_coeff )
           CALL get_aux_coeff(coeffs=local_coeffs,c=c,error=error)
           ub = SIZE ( c, 2 )
           CPPostcondition(ub<=maxgtops,cp_failure_level,routineN,error,failure)
           DO ii = 1, ub
             rep_coeff ( iatom, ii ) = c(iparticle_local,ii )
           END DO
         CASE ( copy_fcoeff )
         CALL get_aux_coeff(coeffs=local_coeffs,fc=fc,error=error)
           CPPostcondition(ub<=maxgtops,cp_failure_level,routineN,error,failure)
           ub = SIZE ( fc, 2 )
           DO ii = 1, ub
             rep_coeff ( iatom, ii ) = fc(iparticle_local, ii )
           END DO
         END SELECT
      END DO
   END DO
 END IF
 CALL mp_sum ( rep_coeff, para_env%group )

 END SUBROUTINE create_replicated_coeff
! *****************************************************************************
!> \brief ...
!> \param rep_coeff ...
!> \param qs_kind_set ...
!> \param local_particles ...
!> \param aux_coeff_set ...
!> \param copy_type_id ...
!> \param copy_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE deallocate_replicated_coeff ( rep_coeff, qs_kind_set, local_particles, &
                                           aux_coeff_set, copy_type_id, copy_out, error )

    REAL(dp), DIMENSION(:, :), POINTER       :: rep_coeff
    TYPE(qs_kind_type), POINTER              :: qs_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    INTEGER, INTENT(IN)                      :: copy_type_id
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_replicated_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ii, ikind, &
                                                iparticle_local, &
                                                nparticle_local, stat
    LOGICAL                                  :: failure, my_copy_out
    REAL(dp), POINTER                        :: c( :, : ), fc( :, : )
    TYPE(aux_coeff_type), POINTER            :: local_coeffs

  NULLIFY ( c, fc )
  my_copy_out = .FALSE.
  failure = .FALSE.
  IF ( PRESENT ( copy_out ) ) my_copy_out = copy_out

  CPPrecondition(ASSOCIATED(rep_coeff),cp_failure_level,routineN,error,failure)
  IF ( .NOT. failure ) THEN
    IF ( my_copy_out ) THEN
      DO ikind = 1, SIZE ( qs_kind_set )
        local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
        IF (.NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
  ! Get the number of particles of ikind ( local )
        nparticle_local = local_particles%n_el(ikind)
  ! Loop over the (local) particles
        DO iparticle_local=1,nparticle_local
           iatom = local_particles%list(ikind)%array(iparticle_local)
  ! Assigning the local coefficients:
           SELECT CASE ( copy_type_id )
           CASE ( copy_coeff )
             CALL get_aux_coeff(coeffs=local_coeffs,c=c,error=error)
             DO ii = 1, SIZE ( rep_coeff, 2 )
               c(iparticle_local,ii) = c(iparticle_local,ii) + rep_coeff ( iatom, ii )
             END DO
           CASE ( copy_fcoeff )
             CALL get_aux_coeff(coeffs=local_coeffs,fc=fc,error=error)
             DO ii = 1, SIZE ( rep_coeff, 2 )
               fc(iparticle_local,ii) = fc(iparticle_local,ii) + rep_coeff ( iatom, ii )
             END DO
           END SELECT
        END DO
      END DO
    END IF
    DEALLOCATE ( rep_coeff, stat=stat )
    CPPostcondition(stat==0,cp_failure_level,routineN,error,failure)
    NULLIFY ( rep_coeff )
  END IF
  END SUBROUTINE deallocate_replicated_coeff
! *****************************************************************************
!> \brief ...
!> \param coeffs ...
!> \param c ...
!> \param fc ...
!> \param n_els ...
!> \param ncoef_atom ...
!> \param ncoef_of_kind ...
!> \param id_nr ...
!> \param ref_count ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_aux_coeff(coeffs, c, fc, n_els, ncoef_atom,ncoef_of_kind,id_nr,ref_count,error)

    TYPE(aux_coeff_type), POINTER            :: coeffs
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: c, fc
    INTEGER, INTENT(OUT), OPTIONAL           :: n_els, ncoef_atom, &
                                                ncoef_of_kind, id_nr, &
                                                ref_count
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_aux_coeff', &
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
       IF (PRESENT(c)) c => coeffs%c
       IF (PRESENT(fc)) fc => coeffs%fc
    END IF

  END SUBROUTINE get_aux_coeff
! *****************************************************************************
!> \brief ...
!> \param aux_coeff_set ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aux_coeff_set_initialize(aux_coeff_set,para_env,error)

    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aux_coeff_set_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, ikind, n_els, &
                                                ncoef_atom, ncoef_atom_max, &
                                                ncoef_of_kind, ncoef_tot, &
                                                nel_tot, nkind
    TYPE(aux_coeff_type), POINTER            :: coeffs

    group=para_env%group
    ! initialize to default values

    nel_tot=0
    ncoef_tot=0
    ncoef_atom_max=0

    ! get total numbers
    nkind=SIZE(aux_coeff_set%coeffs_of_kind)
    DO ikind=1,nkind
       coeffs=> aux_coeff_set%coeffs_of_kind(ikind)%coeffs
       IF(.NOT.ASSOCIATED(coeffs)) CYCLE
       CALL get_aux_coeff(coeffs,n_els=n_els,ncoef_atom=ncoef_atom,&
            ncoef_of_kind=ncoef_of_kind,error=error)
       nel_tot=nel_tot+n_els
       ncoef_tot=ncoef_tot+ncoef_of_kind
       ncoef_atom_max=MAX(ncoef_atom_max,ncoef_atom)
    END DO
    CALL mp_sum(nel_tot,group)
    CALL mp_sum(ncoef_tot,group)
    CALL mp_max(ncoef_atom_max,group)

    aux_coeff_set%nel_tot=nel_tot
    aux_coeff_set%ncoef_tot=ncoef_tot
    aux_coeff_set%ncoef_atom_max=ncoef_atom_max

  END SUBROUTINE aux_coeff_set_initialize

END MODULE scp_coeff_types
