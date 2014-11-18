!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A collection of methods to treat the QM/MM links
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_links_methods
  
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: add_set_type,&
                                             qmmm_env_qm_type,&
                                             qmmm_imomm_link_type,&
                                             qmmm_links_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_links_methods'
  PUBLIC ::  qmmm_link_Imomm_coord,  &
             qmmm_link_Imomm_forces, &
             qmmm_added_chrg_coord,  &
             qmmm_added_chrg_forces

CONTAINS

! *****************************************************************************
!> \brief correct the position for qm/mm IMOMM link type
!> \param qmmm_links ...
!> \param particles ...
!> \param qm_atom_index ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE  qmmm_link_Imomm_coord( qmmm_links, particles, qm_atom_index, error)
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_link_Imomm_coord', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilink, ip, ip_mm, ip_qm, &
                                                mm_index, n_imomm, qm_index
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    TYPE(qmmm_imomm_link_type), POINTER      :: my_link

    failure = .FALSE.
    IF (.NOT.failure) THEN
       n_imomm = SIZE(qmmm_links%imomm)
       CPPrecondition(n_imomm /= 0,cp_failure_level,routineP,error,failure)
       DO ilink = 1, n_imomm
          my_link => qmmm_links%imomm(ilink)%link
          qm_index = my_link%qm_index
          mm_index = my_link%mm_index
          alpha    = 1.0_dp/my_link%alpha
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == qm_index) EXIT
          END DO
          CALL cp_assert(ip/=SIZE(qm_atom_index)+1,cp_failure_level,cp_assertion_failed,routineP,&
               "QM atom index ("//cp_to_string(qm_index)//") specified in the LINK section nr.("//&
               cp_to_string(ilink)//") is not defined as a QM atom! Please inspect your QM_KIND sections. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ip_qm = ip
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == mm_index) EXIT
          END DO
          CALL cp_assert(ip/=SIZE(qm_atom_index)+1,cp_failure_level,cp_assertion_failed,routineP,&
               "Error in setting up the MM atom index ("//cp_to_string(mm_index)//&
               ") specified in the LINK section nr.("//cp_to_string(ilink)//"). Please report this bug! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ip_mm = ip
          particles(ip_mm)%r = alpha *  particles(ip_mm)%r + (1.0_dp-alpha) * particles(ip_qm)%r
       END DO
    END IF

  END SUBROUTINE qmmm_link_Imomm_coord

! *****************************************************************************
!> \brief correct the forces for qm/mm IMOMM link type
!> \param qmmm_links ...
!> \param particles_qm ...
!> \param qm_atom_index ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE  qmmm_link_Imomm_forces(qmmm_links, particles_qm, qm_atom_index, error)
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_link_Imomm_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilink, ip, ip_mm, ip_qm, &
                                                mm_index, n_imomm, qm_index
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    TYPE(qmmm_imomm_link_type), POINTER      :: my_link

    failure = .FALSE.
    IF (.NOT.failure) THEN
       n_imomm = SIZE(qmmm_links%imomm)
       CPPrecondition(n_imomm /= 0,cp_failure_level,routineP,error,failure)
       DO ilink = 1, n_imomm
          my_link => qmmm_links%imomm(ilink)%link
          qm_index = my_link%qm_index
          mm_index = my_link%mm_index
          alpha    = 1.0_dp/my_link%alpha
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == qm_index) EXIT
          END DO
          CALL cp_assert(ip/=SIZE(qm_atom_index)+1,cp_failure_level,cp_assertion_failed,routineP,&
               "QM atom index ("//cp_to_string(qm_index)//") specified in the LINK section nr.("//&
               cp_to_string(ilink)//") is not defined as a QM atom! Please inspect your QM_KIND sections. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ip_qm = ip
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == mm_index) EXIT
          END DO
          CALL cp_assert(ip/=SIZE(qm_atom_index)+1,cp_failure_level,cp_assertion_failed,routineP,&
               "Error in setting up the MM atom index ("//cp_to_string(mm_index)//&
               ") specified in the LINK section nr.("//cp_to_string(ilink)//"). Please report this bug! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ip_mm = ip
          particles_qm(ip_qm)%f = particles_qm(ip_qm)%f + particles_qm(ip_mm)%f * (1.0_dp-alpha)
          particles_qm(ip_mm)%f = particles_qm(ip_mm)%f * alpha
       END DO
    END IF

  END SUBROUTINE qmmm_link_Imomm_forces

! *****************************************************************************
!> \brief correct the position for added charges in qm/mm link scheme
!> \param qmmm_env ...
!> \param particles ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_added_chrg_coord(qmmm_env, particles, error)
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_added_chrg_coord', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, Index1, Index2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    TYPE(add_set_type), POINTER              :: added_charges

    failure = .FALSE.
    added_charges => qmmm_env%added_charges

    DO i = 1, added_charges%num_mm_atoms
       Index1 = added_charges%add_env(i)%Index1
       Index2 = added_charges%add_env(i)%Index1
       alpha  = added_charges%add_env(i)%alpha
       added_charges%added_particles(i)%r = alpha *  particles(Index1)%r + (1.0_dp-alpha) * particles(Index2)%r
    END DO

  END SUBROUTINE qmmm_added_chrg_coord

! *****************************************************************************
!> \brief correct the forces due to the  added charges in qm/mm link scheme
!> \param qmmm_env ...
!> \param particles ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_added_chrg_forces(qmmm_env, particles, error)
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_added_chrg_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, Index1, Index2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    TYPE(add_set_type), POINTER              :: added_charges

    failure = .FALSE.
    added_charges => qmmm_env%added_charges

    DO i = 1, added_charges%num_mm_atoms
       Index1 = added_charges%add_env(i)%Index1
       Index2 = added_charges%add_env(i)%Index1
       alpha  = added_charges%add_env(i)%alpha
       particles(Index1)%f = particles(Index1)%f + alpha           * added_charges%added_particles(i)%f
       particles(Index2)%f = particles(Index2)%f + (1.0_dp -alpha) * added_charges%added_particles(i)%f
    END DO

  END SUBROUTINE qmmm_added_chrg_forces

END MODULE qmmm_links_methods
