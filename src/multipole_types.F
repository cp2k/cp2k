!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Multipole structure: for multipole (fixed and induced) in FF based MD
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
MODULE multipole_types
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE external_potential_types,        ONLY: fist_potential_type,&
                                             get_potential
  USE f77_blas
  USE input_constants,                 ONLY: do_multipole_charge,&
                                             do_multipole_dipole,&
                                             do_multipole_none,&
                                             do_multipole_quadrupole
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: multipole_type,&
            create_multipole_type,&
            release_multipole_type,&
            retain_multipole_type

! *****************************************************************************
!> \brief Define multipole type
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  TYPE multipole_type
     INTEGER :: id_nr, ref_count
     LOGICAL, DIMENSION(3)                    :: task
     REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
     REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: dipoles
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: quadrupoles
  END TYPE multipole_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'multipole_types'
  INTEGER, PRIVATE, SAVE :: last_multipole_id_nr=0

CONTAINS

! *****************************************************************************
!> \brief Create a multipole type
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE create_multipole_type(multipoles, particle_set, subsys_section, ewald_env, error)
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_multipole_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind2, iparticle, j, &
                                                max_multipole, n_rep, &
                                                nparticles, stat
    LOGICAL                                  :: do_multipoles, explicit, &
                                                failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(multipoles),cp_failure_level,routineP,error,failure)
    CALL ewald_env_get(ewald_env, do_multipoles=do_multipoles, max_multipole=max_multipole,&
         error=error)
    IF (do_multipoles) THEN
       ALLOCATE(multipoles, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       last_multipole_id_nr=last_multipole_id_nr+1
       multipoles%id_nr=last_multipole_id_nr
       multipoles%ref_count=1
       multipoles%task = .FALSE.
       NULLIFY(multipoles%charges)
       NULLIFY(multipoles%radii)
       NULLIFY(multipoles%dipoles)
       NULLIFY(multipoles%quadrupoles)
       SELECT CASE(max_multipole)
       CASE(do_multipole_none)
          ! Do nothing..
       CASE(do_multipole_charge)
          multipoles%task(1:1) = .TRUE.
       CASE(do_multipole_dipole)
          multipoles%task(1:2) = .TRUE.
       CASE(do_multipole_quadrupole)
          multipoles%task(1:3) = .TRUE.
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       nparticles = SIZE(particle_set)
       IF (multipoles%task(1)) THEN
          ALLOCATE(multipoles%charges(nparticles),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(multipoles%radii(nparticles),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Fill in charge array
          DO iparticle = 1, nparticles
             !atomic_kind =>
             CALL get_atomic_kind(particle_set(iparticle)%atomic_kind,&
                                  fist_potential=fist_potential)
             CALL get_potential(fist_potential,qeff=multipoles%charges(iparticle),&
                                mm_radius=multipoles%radii(iparticle))
          END DO
       END IF
       IF (multipoles%task(2)) THEN
          ALLOCATE(multipoles%dipoles(3,nparticles),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Fill in dipole array (if specified)
          work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%DIPOLES",error=error)
          CALL section_vals_get(work_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",n_rep_val=n_rep, error=error)
             CPPostcondition(n_rep==nparticles,cp_failure_level,routineP,error,failure)
             DO iparticle = 1, n_rep
                CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",i_rep_val=iparticle, r_vals=work, error=error)
                multipoles%dipoles(1:3,iparticle) = work
             END DO
          ELSE
             multipoles%dipoles = 0.0_dp
          END IF
       END IF
       IF (multipoles%task(3)) THEN
          ALLOCATE(multipoles%quadrupoles(3,3,nparticles),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Fill in quadrupole array (if specified)
          work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%QUADRUPOLES",error=error)
          CALL section_vals_get(work_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",n_rep_val=n_rep, error=error)
             CPPostcondition(n_rep==nparticles,cp_failure_level,routineP,error,failure)
             DO iparticle = 1, n_rep
                CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",i_rep_val=iparticle, r_vals=work, error=error)
                DO i = 1, 3
                   DO j = 1, 3
                      ind2 = 3*(MIN(i,j)-1)-(MIN(i,j)*(MIN(i,j)-1))/2+MAX(i,j)
                      multipoles%quadrupoles(i,j,iparticle) = work(ind2)
                   END DO
                END DO
             END DO
          ELSE
             multipoles%quadrupoles = 0.0_dp
          END IF
       END IF
    END IF
  END SUBROUTINE create_multipole_type

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE release_multipole_type(multipoles, error)
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_multipole_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure =.FALSE.
    IF (ASSOCIATED(multipoles)) THEN
       CPPostcondition(multipoles%ref_count>0,cp_failure_level,routineP,error,failure)
       multipoles%ref_count=multipoles%ref_count-1
       IF (multipoles%ref_count==0) THEN
          IF (ASSOCIATED(multipoles%charges)) THEN
             DEALLOCATE(multipoles%charges, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(multipoles%radii)) THEN
             DEALLOCATE(multipoles%radii, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(multipoles%dipoles)) THEN
             DEALLOCATE(multipoles%dipoles, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(multipoles%quadrupoles)) THEN
             DEALLOCATE(multipoles%quadrupoles, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(multipoles, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE release_multipole_type

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE retain_multipole_type(multipoles, error)
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_multipole_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(multipoles)) THEN
       CPPreconditionNoFail(multipoles%ref_count>0,cp_failure_level,routineP,error)
       multipoles%ref_count=multipoles%ref_count+1
    END IF
  END SUBROUTINE retain_multipole_type

END MODULE multipole_types
