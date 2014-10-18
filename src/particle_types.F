!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the data structure for the particle information.
!> \par History
!>      - Atomic kind added in particle_type (MK,08.01.2002)
!>      - Functionality for particle_type added (MK,14.01.2002)
!>      - Allow for general coordinate input (MK,13.09.2003)
!>      - Molecule concept introduced (MK,26.09.2003)
!>      - Last atom information added (jgh,23.05.2004)
!>      - particle_type cleaned (MK,03.02.2005)
!> \author CJM, MK
! *****************************************************************************
MODULE particle_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'particle_types'

  ! Data types
! *****************************************************************************
  TYPE particle_type
     TYPE(atomic_kind_type), POINTER       :: atomic_kind ! atomic kind information
     REAL(KIND = dp), DIMENSION(3)         :: f,&         ! force
                                              r,&         ! position
                                              v           ! velocity
     ! Particle dependent terms for shell-model
     INTEGER                               :: atom_index, t_region_index, shell_index
  END TYPE particle_type

  ! Public data types

  PUBLIC :: particle_type

  ! Public subroutines

  PUBLIC :: allocate_particle_set,&
            deallocate_particle_set,&
            get_particle_set,&
            update_particle_set,&
            update_particle_pos_or_vel,&
            get_particle_pos_or_vel

CONTAINS

! *****************************************************************************
!> \brief   Allocate a particle set.
!> \param particle_set ...
!> \param nparticle ...
!> \param error ...
!> \date    14.01.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_particle_set(particle_set,nparticle,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: nparticle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(particle_set)) THEN
       CALL deallocate_particle_set(particle_set,error=error)
    END IF
    ALLOCATE (particle_set(nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iparticle=1,nparticle
      NULLIFY (particle_set(iparticle)%atomic_kind)
      particle_set(iparticle)%f(:) = 0.0_dp
      particle_set(iparticle)%r(:) = 0.0_dp
      particle_set(iparticle)%v(:) = 0.0_dp
      particle_set(iparticle)%shell_index = 0
      particle_set(iparticle)%atom_index = 0
      particle_set(iparticle)%t_region_index = 0
    END DO

  END SUBROUTINE allocate_particle_set

! *****************************************************************************
!> \brief   Deallocate a particle set.
!> \param particle_set ...
!> \param error ...
!> \date    14.01.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_particle_set(particle_set,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(particle_set)) THEN
      DEALLOCATE (particle_set,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer particle_set is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_particle_set

! *****************************************************************************
!> \brief   Get the components of a particle set.
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param first_sgf ...
!> \param last_sgf ...
!> \param nsgf ...
!> \param basis_set_id ...
!> \param basis ...
!> \param error ...
!> \date    14.01.2002
!> \par History
!>      - particle type cleaned (13.10.2003,MK)
!>      - refactoring and add besis set option (17.08.2010,jhu)
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_particle_set(particle_set,qs_kind_set,first_sgf,last_sgf,nsgf,&
                              basis_set_id,basis,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER, DIMENSION(:), INTENT(INOUT), &
      OPTIONAL                               :: first_sgf, last_sgf, nsgf
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, iparticle, isgf, &
                                                nparticle, ns
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)

    nparticle = SIZE(particle_set)
    IF (PRESENT(first_sgf)) THEN
      CPPrecondition(SIZE(first_sgf) >= nparticle,cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(last_sgf)) THEN
      CPPrecondition(SIZE(last_sgf) >= nparticle,cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(nsgf)) THEN
      CPPrecondition(SIZE(nsgf) >= nparticle,cp_failure_level,routineP,error,failure)
    END IF

    IF (PRESENT(first_sgf).OR.PRESENT(last_sgf).OR.PRESENT(nsgf)) THEN
      isgf = 0
      DO iparticle=1,nparticle
        CALL get_atomic_kind(particle_set(iparticle)%atomic_kind, kind_number=ikind)
        IF (PRESENT(basis)) THEN
          CALL get_gto_basis_set(gto_basis_set=basis(ikind)%gto_basis_set,nsgf=ns)
        ELSEIF (PRESENT(basis_set_id)) THEN
          CALL get_qs_kind(qs_kind_set(ikind), nsgf=ns, basis_set_id=basis_set_id)
        ELSE
          CALL get_qs_kind(qs_kind_set(ikind), nsgf=ns)
        END IF
        IF (PRESENT(nsgf)) nsgf(iparticle) = ns
        IF (PRESENT(first_sgf)) first_sgf(iparticle) = isgf + 1
        isgf = isgf + ns
        IF (PRESENT(last_sgf)) last_sgf(iparticle) = isgf
      END DO
    END IF
    IF (PRESENT(first_sgf)) THEN
      IF (SIZE(first_sgf) > nparticle) first_sgf(nparticle+1) = isgf + 1
    END IF

  END SUBROUTINE get_particle_set


! *****************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param int_group ...
!> \param pos ...
!> \param vel ...
!> \param for ...
!> \param add ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE update_particle_set ( particle_set, int_group,  pos, vel, for, add, error )

    TYPE(particle_type), POINTER             :: particle_set( : )
    INTEGER, INTENT(IN)                      :: int_group
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: pos( :, : ), vel( :, : ), &
                                                for(:,:)
    LOGICAL, INTENT(IN), OPTIONAL            :: add
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iparticle, nparticle
    LOGICAL                                  :: failure, my_add, update_for, &
                                                update_pos, update_vel

    CALL timeset(routineN,handle)

    failure = .FALSE.
    nparticle  = SIZE ( particle_set )
    update_pos = PRESENT(pos)
    update_vel = PRESENT(vel)
    update_for = PRESENT(for)
    my_add     = .FALSE.
    IF (PRESENT(add)) my_add = add

    IF (update_pos) THEN
       CALL mp_sum( pos,int_group)
       IF (my_add) THEN
          DO iparticle=1,nparticle
             particle_set(iparticle)%r(:) = particle_set(iparticle)%r(:) + pos(:,iparticle)
          END DO
       ELSE
          DO iparticle=1,nparticle
             particle_set(iparticle)%r(:) = pos(:,iparticle)
          END DO
       END IF
    END IF
    IF (update_vel) THEN
       CALL mp_sum( vel,int_group)
       IF (my_add) THEN
          DO iparticle=1,nparticle
             particle_set(iparticle)%v(:) = particle_set(iparticle)%v(:) + vel(:,iparticle)
          END DO
       ELSE
          DO iparticle=1,nparticle
             particle_set(iparticle)%v(:) = vel(:,iparticle)
          END DO
       END IF
    END IF
    IF (update_for) THEN
       CALL mp_sum( for,int_group)
       IF (my_add) THEN
          DO iparticle=1,nparticle
             particle_set(iparticle)%f(:) = particle_set(iparticle)%f(:) + for(:,iparticle)
          END DO
       ELSE
          DO iparticle=1,nparticle
             particle_set(iparticle)%f(:) = for(:,iparticle)
          END DO
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE update_particle_set

! *****************************************************************************
!> \brief   Return the atomic position or velocity of atom iatom in x from a
!>          packed vector even if core-shell particles are present
!> \param iatom ...
!> \param particle_set ...
!> \param vector ...
!> \retval x ...
!> \date    25.11.2010
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION get_particle_pos_or_vel(iatom,particle_set,vector) RESULT(x)

    INTEGER, INTENT(IN)                      :: iatom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector
    REAL(KIND=dp), DIMENSION(3)              :: x

    INTEGER                                  :: ic, is
    REAL(KIND=dp)                            :: fc, fs, mass

    ic = 3*(iatom - 1)
    IF (particle_set(iatom)%shell_index == 0) THEN
       x(1:3) = vector(ic+1:ic+3)
    ELSE
       is = 3*(SIZE(particle_set) + particle_set(iatom)%shell_index - 1)
       mass = particle_set(iatom)%atomic_kind%mass
       fc = particle_set(iatom)%atomic_kind%shell%mass_core/mass
       fs = particle_set(iatom)%atomic_kind%shell%mass_shell/mass
       x(1:3) = fc*vector(ic+1:ic+3) + fs*vector(is+1:is+3)
    END IF

  END FUNCTION get_particle_pos_or_vel

! *****************************************************************************
!> \brief   Update the atomic position or velocity by x and return the updated
!>          atomic position or velocity in x even if core-shell particles are
!>          present
!> \param iatom ...
!> \param particle_set ...
!> \param x ...
!> \param vector ...
!> \date    26.11.2010
!> \author  Matthias Krack
!> \version 1.0
!> \note    particle-set is not changed, only the positions or velocities in
!>          the packed vector are updated
! *****************************************************************************
  SUBROUTINE update_particle_pos_or_vel(iatom,particle_set,x,vector)

    INTEGER, INTENT(IN)                      :: iatom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: x
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: vector

    INTEGER                                  :: ic, is
    REAL(KIND=dp)                            :: fc, fs, mass

    ic = 3*(iatom - 1)
    IF (particle_set(iatom)%shell_index == 0) THEN
       vector(ic+1:ic+3) = vector(ic+1:ic+3) + x(1:3)
       x(1:3) = vector(ic+1:ic+3)
    ELSE
       is = 3*(SIZE(particle_set) + particle_set(iatom)%shell_index - 1)
       mass = particle_set(iatom)%atomic_kind%mass
       fc = particle_set(iatom)%atomic_kind%shell%mass_core/mass
       fs = particle_set(iatom)%atomic_kind%shell%mass_shell/mass
       vector(ic+1:ic+3) = vector(ic+1:ic+3) + x(1:3)
       vector(is+1:is+3) = vector(is+1:is+3) + x(1:3)
       x(1:3) = fc*vector(ic+1:ic+3) + fs*vector(is+1:is+3)
    END IF

  END SUBROUTINE update_particle_pos_or_vel

END MODULE particle_types
