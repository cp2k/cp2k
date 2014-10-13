!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Data type and methods dealing with PI calcs in normal mode coords
!> \author fawzi
!> \par    History
!>         2006-02 created
!>         2006-11 modified so it might actually work [hforbert]
!>         2009-04-07 moved from pint_types module to a separate file [lwalewski]
! *****************************************************************************
MODULE pint_normalmode
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE pint_types,                      ONLY: normalmode_env_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_normalmode'

  INTEGER, SAVE, PRIVATE :: last_normalmode_id=0

  PUBLIC :: normalmode_env_create
  PUBLIC :: normalmode_release
  PUBLIC :: normalmode_init_masses
  PUBLIC :: normalmode_x2u
  PUBLIC :: normalmode_u2x
  PUBLIC :: normalmode_f2uf
  PUBLIC :: normalmode_calc_uf_h

CONTAINS

! ***************************************************************************
!> \brief creates the data needed for a normal mode transformation
!> \param normalmode_env ...
!> \param normalmode_section ...
!> \param p ...
!> \param kT ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_env_create(normalmode_env,normalmode_section,p,kT,error)
    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    TYPE(section_vals_type), POINTER         :: normalmode_section
    INTEGER, INTENT(in)                      :: p
    REAL(kind=dp), INTENT(in)                :: kT
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, li, stat
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(p)              :: lambda

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(normalmode_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
      ALLOCATE(normalmode_env,stat=stat)
      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
      ALLOCATE(normalmode_env%x2u(p,p),stat=stat)
      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
      ALLOCATE(normalmode_env%u2x(p,p),stat=stat)
      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
     last_normalmode_id=last_normalmode_id+1
     normalmode_env%id_nr=last_normalmode_id
     normalmode_env%ref_count=1

     normalmode_env%p=p

     CALL section_vals_val_get(normalmode_section,"Q_CENTROID", &
                               r_val=normalmode_env%Q_centroid,error=error)
     CALL section_vals_val_get(normalmode_section,"Q_BEAD", &
                               r_val=normalmode_env%Q_bead,error=error)
     CALL section_vals_val_get(normalmode_section,"MODEFACTOR", &
                               r_val=normalmode_env%modefactor,error=error)

     IF (normalmode_env%Q_centroid < 0.0_dp) THEN
        normalmode_env%Q_centroid = -normalmode_env%Q_centroid/(kT*p)
     END IF
     IF (normalmode_env%Q_bead < 0.0_dp) THEN
        normalmode_env%Q_bead = -normalmode_env%Q_bead/(kT*p)
     END IF

     normalmode_env%harm = p*kT*kT/normalmode_env%modefactor

     ! set up the transformation matrices
     DO i = 1, p
        lambda(i) = 2.0_dp*(1.0_dp - COS(pi*(i/2)*2.0_dp/p))
        DO j = 1, p
           k = ((i/2)*(j-1))/p
           k = (i/2)*(j-1)-k*p
           li = 2*(i-2*(i/2))*p-p
           normalmode_env%u2x(j,i)=SQRT(2.0_dp/p)*SIN(twopi*(k+0.125_dp*li)/p)
        END DO
     END DO
     lambda(1) = 1.0_dp/(p*normalmode_env%modefactor)
     DO i = 1, p
        DO j = 1, p
           normalmode_env%x2u(i,j)=SQRT(lambda(i)*normalmode_env%modefactor)*&
                                   normalmode_env%u2x(j,i)
        END DO
     END DO
     DO i = 1, p
        DO j = 1, p
           normalmode_env%u2x(i,j)=normalmode_env%u2x(i,j)/ &
                                   SQRT(lambda(j)*normalmode_env%modefactor)
        END DO
     END DO
    END IF
  END SUBROUTINE normalmode_env_create

! ***************************************************************************
!> \brief releases the normalmode environment
!> \param normalmode_env the normalmode_env to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_release(normalmode_env,error)

    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(normalmode_env)) THEN
      CPPrecondition(normalmode_env%ref_count>0,cp_failure_level,routineP,error,failure)
      normalmode_env%ref_count=normalmode_env%ref_count-1
      IF (normalmode_env%ref_count==0) THEN
        DEALLOCATE(normalmode_env%x2u,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(normalmode_env%u2x,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(normalmode_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
    END IF
    NULLIFY(normalmode_env)

    RETURN
  END SUBROUTINE normalmode_release

! ***************************************************************************
!> \brief retains a normalmode_env
!> \param normalmode_env the normalmode_env to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_retain(normalmode_env,error)

    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(normalmode_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CPPrecondition(normalmode_env%ref_count>0,cp_failure_level,routineP,error,failure)
      normalmode_env%ref_count=normalmode_env%ref_count+1
    END IF

    RETURN
  END SUBROUTINE normalmode_retain

! ***************************************************************************
!> \brief initializes the masses and fictitious masses compatible with the
!>      normal mode information
!> \param normalmode_env the definition of the normal mode transformation
!> \param mass *input* the masses of the particles
!> \param mass_beads masses of the beads
!> \param mass_fict the fictitious masses
!> \param Q masses of the nose thermostats
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_init_masses(normalmode_env,mass,mass_beads,mass_fict,&
     Q,error)

    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    REAL(kind=dp), DIMENSION(:), INTENT(in)  :: mass
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL                  :: mass_beads, mass_fict
    REAL(kind=dp), DIMENSION(:), &
      INTENT(out), OPTIONAL                  :: Q
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_init_masses', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, ib
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT. failure) THEN
      IF (PRESENT(Q)) THEN
        Q=normalmode_env%Q_bead
        Q(1)=normalmode_env%Q_centroid
      END IF
      IF (PRESENT(mass_beads).OR.PRESENT(mass_fict)) THEN
         IF (PRESENT(mass_beads)) THEN
            DO iat=1,SIZE(mass)
               mass_beads(1,iat)=0.0_dp
               DO ib=2,normalmode_env%p
                  mass_beads(ib,iat)=mass(iat)
               END DO
            END DO
         END IF
         IF (PRESENT(mass_fict)) THEN
            DO iat=1,SIZE(mass)
               DO ib=1,normalmode_env%p
                  mass_fict(ib,iat)=mass(iat)
               END DO
            END DO
         END IF
      END IF
    END IF

    RETURN
  END SUBROUTINE normalmode_init_masses

! ***************************************************************************
!> \brief Transforms from the x into the u variables using a normal mode
!>      transformation for the positions
!> \param normalmode_env the environment for the normal mode transformation
!> \param ux will contain the u variable
!> \param x the positions to transform
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_x2u(normalmode_env,ux,x,error)
    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_x2u', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(normalmode_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(normalmode_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL DGEMM('N','N',normalmode_env%p,SIZE(x,2),normalmode_env%p,1.0_dp,&
        normalmode_env%x2u(1,1),SIZE(normalmode_env%x2u,1),x(1,1),SIZE(x,1),&
        0.0_dp,ux,SIZE(ux,1))
    END IF
    RETURN
  END SUBROUTINE normalmode_x2u

! ***************************************************************************
!> \brief transform from the u variable to the x (back normal mode
!>      transformation for the positions)
!> \param normalmode_env the environment for the normal mode transformation
!> \param ux the u variable (positions to be backtransformed)
!> \param x will contain the positions
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_u2x(normalmode_env,ux,x,error)
    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_u2x', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(normalmode_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(normalmode_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL DGEMM('N','N',normalmode_env%p,SIZE(ux,2),normalmode_env%p,1.0_dp,&
        normalmode_env%u2x(1,1),SIZE(normalmode_env%u2x,1),ux(1,1),SIZE(ux,1),&
        0.0_dp,x,SIZE(x,1))
    END IF
    RETURN
  END SUBROUTINE normalmode_u2x

! ***************************************************************************
!> \brief normalmode transformation for the forces
!> \param normalmode_env the environment for the normal mode transformation
!> \param uf will contain the forces for the transformed variables afterwards
!> \param f the forces to transform
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_f2uf(normalmode_env,uf,f,error)
    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: uf
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: f
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_f2uf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(normalmode_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(normalmode_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL DGEMM('T','N',normalmode_env%p,SIZE(f,2),normalmode_env%p,1.0_dp,&
        normalmode_env%u2x(1,1),SIZE(normalmode_env%u2x,1),f(1,1),SIZE(f,1),&
        0.0_dp,uf,SIZE(uf,1))
    END IF
    RETURN
  END SUBROUTINE normalmode_f2uf

! ***************************************************************************
!> \brief calculates the harmonic force in the normal mode basis
!> \param normalmode_env the normal mode environment
!> \param mass_beads the masses of the beads
!> \param ux the positions of the beads in the staging basis
!> \param uf_h the harmonic forces (not accelerations)
!> \param e_h ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Harald Forbert
! *****************************************************************************
  SUBROUTINE normalmode_calc_uf_h(normalmode_env,mass_beads,ux,uf_h,e_h,error)
    TYPE(normalmode_env_type), POINTER       :: normalmode_env
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: mass_beads, ux, uf_h
    REAL(KIND=dp), INTENT(OUT)               :: e_h
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'normalmode_calc_uf_h', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ibead, idim
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: f

    failure=.FALSE.

    IF (.NOT. failure) THEN
      e_h = 0.0_dp
      DO idim = 1, SIZE(mass_beads,2)

        ! starting at 2 since the centroid is at 1 and it's mass_beads
        ! SHOULD be zero anyways:

        uf_h(1,idim) = 0.0_dp
        DO ibead = 2, normalmode_env%p
           f = - mass_beads(ibead,idim)*normalmode_env%harm*ux(ibead,idim)
           uf_h(ibead,idim)= f
           ! - to cancel the - in the force f.
           e_h = e_h - 0.5_dp*ux(ibead,idim)*f
        END DO

      END DO
    END IF
    RETURN
  END SUBROUTINE normalmode_calc_uf_h

END MODULE pint_normalmode
