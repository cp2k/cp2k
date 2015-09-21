!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Multipole structure: for multipole (fixed and induced) in FF based MD
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! *****************************************************************************
MODULE multipole_types
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE external_potential_types,        ONLY: fist_potential_type,&
                                             get_potential
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: multipole_type,&
            create_multipole_type,&
            release_multipole_type,&
            retain_multipole_type

  INTEGER, PARAMETER, PUBLIC               :: do_multipole_none=-1,&
                                              do_multipole_charge=0,&
                                              do_multipole_dipole=1,&
                                              do_multipole_quadrupole=2

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
!> \param multipoles ...
!> \param particle_set ...
!> \param subsys_section ...
!> \param max_multipole ...
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE create_multipole_type(multipoles, particle_set, subsys_section, max_multipole)
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    INTEGER, INTENT(IN)                      :: max_multipole

    CHARACTER(len=*), PARAMETER :: routineN = 'create_multipole_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind2, iparticle, j, n_rep, &
                                                nparticles
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    ALLOCATE(multipoles)

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
       CPASSERT(.FALSE.)
    END SELECT
    nparticles = SIZE(particle_set)
    IF (multipoles%task(1)) THEN
       ALLOCATE(multipoles%charges(nparticles))
       ALLOCATE(multipoles%radii(nparticles))
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
       ALLOCATE(multipoles%dipoles(3,nparticles))
       ! Fill in dipole array (if specified)
       work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%DIPOLES")
       CALL section_vals_get(work_section, explicit=explicit)
       IF (explicit) THEN
          CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",n_rep_val=n_rep)
          CPASSERT(n_rep==nparticles)
          DO iparticle = 1, n_rep
             CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",i_rep_val=iparticle, r_vals=work)
             multipoles%dipoles(1:3,iparticle) = work
          END DO
       ELSE
          multipoles%dipoles = 0.0_dp
       END IF
    END IF
    IF (multipoles%task(3)) THEN
       ALLOCATE(multipoles%quadrupoles(3,3,nparticles))
       ! Fill in quadrupole array (if specified)
       work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%QUADRUPOLES")
       CALL section_vals_get(work_section, explicit=explicit)
       IF (explicit) THEN
          CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",n_rep_val=n_rep)
          CPASSERT(n_rep==nparticles)
          DO iparticle = 1, n_rep
             CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",i_rep_val=iparticle, r_vals=work)
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
  END SUBROUTINE create_multipole_type

! *****************************************************************************
!> \brief ...
!> \param multipoles ...
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE release_multipole_type(multipoles)
    TYPE(multipole_type), POINTER            :: multipoles

    CHARACTER(len=*), PARAMETER :: routineN = 'release_multipole_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.
    IF (ASSOCIATED(multipoles)) THEN
       CPASSERT(multipoles%ref_count>0)
       multipoles%ref_count=multipoles%ref_count-1
       IF (multipoles%ref_count==0) THEN
          IF (ASSOCIATED(multipoles%charges)) THEN
             DEALLOCATE(multipoles%charges)
          END IF
          IF (ASSOCIATED(multipoles%radii)) THEN
             DEALLOCATE(multipoles%radii)
          END IF
          IF (ASSOCIATED(multipoles%dipoles)) THEN
             DEALLOCATE(multipoles%dipoles)
          END IF
          IF (ASSOCIATED(multipoles%quadrupoles)) THEN
             DEALLOCATE(multipoles%quadrupoles)
          END IF
          DEALLOCATE(multipoles)
       END IF
    END IF
  END SUBROUTINE release_multipole_type

! *****************************************************************************
!> \brief ...
!> \param multipoles ...
!> \par History
!>      12.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE retain_multipole_type(multipoles)
    TYPE(multipole_type), POINTER            :: multipoles

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_multipole_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(multipoles)) THEN
       CPASSERT(multipoles%ref_count>0)
       multipoles%ref_count=multipoles%ref_count+1
    END IF
  END SUBROUTINE retain_multipole_type

END MODULE multipole_types
