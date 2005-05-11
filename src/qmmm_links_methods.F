!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_links_methods [1.0] *
!!
!!   NAME MODULE
!!     qmmm_links_methods
!!
!!   FUNCTION
!!     A collection of methods to treat the QM/MM links
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_links_methods
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_links_type,&
                                             qmmm_imomm_link_type,&
                                             qmmm_pseudo_link_type,&
                                             qmmm_env_qm_type,&
                                             add_set_type

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_link_methods'
  PUBLIC ::  qmmm_link_Imomm_coord,  & 
             qmmm_link_Imomm_forces, &
             qmmm_added_chrg_coord,  &
             qmmm_added_chrg_forces


CONTAINS
!!****f* qmmm_links_methods/ qmmm_link_Imomm_coord [1.0] *
!!
!!   NAME
!!     qmmm_link_Imomm_coord
!!
!!   FUNCTION
!!     correct the position for qm/mm IMOMM link type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_link_Imomm_coord( qmmm_links, particles, qm_atom_index, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(particle_type), DIMENSION(:), POINTER   :: particles
    TYPE(qmmm_links_type), POINTER               :: qmmm_links
    INTEGER, DIMENSION(:), POINTER               :: qm_atom_index
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='qmmm_link_Imomm_coord',&
       routineP=moduleN//':'//routineN
    TYPE(qmmm_imomm_link_type), POINTER :: my_link
    LOGICAL :: failure
    INTEGER :: n_imomm, qm_index, mm_index, ilink, ip, ip_qm, ip_mm
    REAL(KIND=dp) :: alpha

    failure = .false.
    IF (.NOT.failure) THEN
       n_imomm = SIZE(qmmm_links%imomm)
       CPPrecondition(n_imomm /= 0,cp_failure_level,routineP,error,failure)
       DO ilink = 1, n_imomm
          my_link => qmmm_links%imomm(ilink)%link
          qm_index = my_link%qm_index
          mm_index = my_link%mm_index
          alpha    = 1.0_dp/my_link%alpha
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == qm_index) exit
          END DO
          CPPrecondition(ip /= SIZE(qm_atom_index)+1,cp_failure_level,routineP,error,failure)
          ip_qm = ip
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == mm_index) exit
          END DO
          CPPrecondition(ip /= SIZE(qm_atom_index)+1,cp_failure_level,routineP,error,failure)
          ip_mm = ip
          particles(ip_mm)%r = alpha *  particles(ip_mm)%r + (1.0_dp-alpha) * particles(ip_qm)%r 
       END DO
    END IF

  END SUBROUTINE qmmm_link_Imomm_coord

!!****f* qmmm_links_methods/ qmmm_link_Imomm_forces [1.0] *
!!
!!   NAME
!!     qmmm_link_Imomm_forces
!!
!!   FUNCTION
!!     correct the forces for qm/mm IMOMM link type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_link_Imomm_forces(qmmm_links, particles_qm, qm_atom_index, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(particle_type), DIMENSION(:), POINTER   :: particles_qm
    TYPE(qmmm_links_type), POINTER               :: qmmm_links
    INTEGER,  DIMENSION(:), POINTER              :: qm_atom_index
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='qmmm_link_Imomm_forces',&
       routineP=moduleN//':'//routineN
    TYPE(qmmm_imomm_link_type), POINTER :: my_link
    LOGICAL :: failure
    INTEGER :: n_imomm, qm_index, mm_index, ilink, ip, ip_qm, ip_mm
    REAL(KIND=dp) :: alpha

    failure = .false.
    IF (.NOT.failure) THEN
       n_imomm = SIZE(qmmm_links%imomm)
       CPPrecondition(n_imomm /= 0,cp_failure_level,routineP,error,failure)
       DO ilink = 1, n_imomm
          my_link => qmmm_links%imomm(ilink)%link
          qm_index = my_link%qm_index
          mm_index = my_link%mm_index
          alpha    = 1.0_dp/my_link%alpha
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == qm_index) exit
          END DO
          CPPrecondition(ip /= SIZE(qm_atom_index)+1,cp_failure_level,routineP,error,failure)
          ip_qm = ip
          DO ip = 1, SIZE(qm_atom_index)
             IF (qm_atom_index(ip) == mm_index) exit
          END DO
          CPPrecondition(ip /= SIZE(qm_atom_index)+1,cp_failure_level,routineP,error,failure)
          ip_mm = ip
          particles_qm(ip_qm)%f = particles_qm(ip_qm)%f + particles_qm(ip_mm)%f * (1.0_dp-alpha)
          particles_qm(ip_mm)%f = particles_qm(ip_mm)%f * alpha          
       END DO
    END IF

  END SUBROUTINE qmmm_link_Imomm_forces

!!****f* qmmm_links_methods/ qmmm_added_chrg_coord [1.0] *
!!
!!   NAME
!!     qmmm_added_chrg_coord
!!
!!   FUNCTION
!!     correct the position for added charges in qm/mm link scheme
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [tlaino]
!!
!!*** **********************************************************************    
  SUBROUTINE qmmm_added_chrg_coord(qmmm_env, particles, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), POINTER   :: particles
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='qmmm_added_chrg_coord',&
       routineP=moduleN//':'//routineN
    LOGICAL :: failure
    TYPE(add_set_type), POINTER                             :: added_charges
    INTEGER :: Index1, Index2,I
    REAL(KIND=dp) :: alpha

    failure = .FALSE.
    added_charges => qmmm_env%added_charges

    DO i = 1, added_charges%num_mm_atoms
       Index1 = added_charges%add_env(i)%Index1
       Index2 = added_charges%add_env(i)%Index1
       alpha  = added_charges%add_env(i)%alpha
       added_charges%added_particles(i)%r = alpha *  particles(Index1)%r + (1.0_dp-alpha) * particles(Index2)%r
    END DO

  END SUBROUTINE qmmm_added_chrg_coord

!!****f* qmmm_links_methods/ qmmm_added_chrg_forces [1.0] *
!!
!!   NAME
!!     qmmm_added_chrg_forces
!!
!!   FUNCTION
!!     correct the forces due to the  added charges in qm/mm link scheme
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [tlaino]
!!
!!*** **********************************************************************    
  SUBROUTINE qmmm_added_chrg_forces(qmmm_env, particles, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), POINTER   :: particles
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN='qmmm_added_chrg_coord',&
       routineP=moduleN//':'//routineN
    LOGICAL :: failure
    TYPE(add_set_type), POINTER                             :: added_charges
    INTEGER :: Index1, Index2,I
    REAL(KIND=dp) :: alpha

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
