!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods to include the effect of an external potential during an MD
!>        or energy calculation
!> \author Teodoro Laino (03.2008) [tlaino]
! *****************************************************************************
MODULE external_potential_methods
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set,&
                                             force_env_type
  USE force_fields_util,               ONLY: get_generic_info
  USE fparser,                         ONLY: evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_list_types,             ONLY: particle_list_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_methods'
  PUBLIC :: add_external_potential

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param error ...
!> \date 03.2008
!> \author Teodoro Laino - University of Zurich [tlaino]
! *****************************************************************************
  SUBROUTINE add_external_potential(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_external_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par
    INTEGER                                  :: a_var, handle, i, iatom, j, &
                                                k, n_var, natom, rep, stat
    INTEGER, DIMENSION(:), POINTER           :: iatms, nparticle
    LOGICAL                                  :: failure, useall
    REAL(KIND=dp)                            :: dedf, dx, energy, err, lerr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ext_pot_section

    failure = .FALSE.
    useall  = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(my_par, my_val, logger, subsys, particles, ext_pot_section, nparticle)
    ext_pot_section => section_vals_get_subs_vals(force_env%force_env_section,&
         "EXTERNAL_POTENTIAL",error=error)
    CALL section_vals_get(ext_pot_section,n_repetition=n_var,error=error)
    DO rep=1, n_var
       natom       =  0
       logger => cp_error_get_logger(error)
       CALL section_vals_val_get(ext_pot_section,"DX",r_val=dx,i_rep_section=rep,error=error)
       CALL section_vals_val_get(ext_pot_section,"ERROR_LIMIT",r_val=lerr,i_rep_section=rep,error=error)
       CALL get_generic_info(ext_pot_section, "FUNCTION", coupling_function, my_par, my_val,&
            input_variables=(/"X","Y","Z"/), i_rep_sec=rep,error=error)
       CALL initf(1)
       CALL parsef(1,TRIM(coupling_function),my_par)

       ! Apply potential on all atoms, computing energy and forces
       NULLIFY(particles, subsys)
       CALL force_env_get(force_env, subsys=subsys, error=error)
       CALL cp_subsys_get(subsys, particles=particles, error=error)
       CALL force_env_get(force_env, additional_potential=energy, error=error)
       CALL section_vals_val_get(ext_pot_section,"ATOMS_LIST",n_rep_val=a_var,i_rep_section=rep,error=error)
       DO k = 1, a_var
          CALL section_vals_val_get(ext_pot_section,"ATOMS_LIST",i_rep_val=k,i_vals=iatms,i_rep_section=rep,&
               error=error)
          CALL reallocate(nparticle,1, natom+SIZE(iatms))
          nparticle(natom+1:natom+SIZE(iatms)) = iatms
          natom = natom + SIZE(iatms)
       END DO
       IF (a_var==0) THEN
          natom  = particles%n_els
          useall = .TRUE.
       END IF
       DO i = 1, natom
          IF (useall) THEN
             iatom = i
          ELSE
             iatom = nparticle(i)
          END IF
          my_val(1)=particles%els(iatom)%r(1)
          my_val(2)=particles%els(iatom)%r(2)
          my_val(3)=particles%els(iatom)%r(3)

          energy = energy + evalf(1,my_val)
          DO j = 1, 3
             dedf = evalfd(1,j,my_val,dx,err)
             IF (ABS(err)>lerr) THEN
                WRITE(this_error,"(A,G12.6,A)")"(",err,")"
                WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
                CALL compress(this_error,.TRUE.)
                CALL compress(def_error,.TRUE.)
                CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
                     'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                     ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                     TRIM(def_error)//' .',error=error,only_ionode=.TRUE.)
             END IF
             particles%els(iatom)%f(j)=particles%els(iatom)%f(j)-dedf
          END DO
       END DO
       CALL force_env_set(force_env, additional_potential=energy, error=error)
       DEALLOCATE(my_par,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_val,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (a_var/=0) THEN
          DEALLOCATE (nparticle,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL finalizef()
    END DO
    CALL timestop(handle)
  END SUBROUTINE add_external_potential

END MODULE external_potential_methods
