!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  USE particle_list_types,             ONLY: particle_list_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_methods'
  PUBLIC :: add_external_potential

CONTAINS

! *****************************************************************************
!> \brief 
!> \author Teodoro Laino - University of Zurich [tlaino] 
!> \date 03.2008
! *****************************************************************************
  SUBROUTINE add_external_potential(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_external_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par, par
    INTEGER                                  :: handle, iatom, j, n_rep, &
                                                output_unit, stat, sz_par
    LOGICAL                                  :: explicit, failure, ionode
    REAL(KIND=dp)                            :: dedf, dx, energy, err, lerr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val, val
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ext_pot_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    ext_pot_section => section_vals_get_subs_vals(force_env%force_env_section,&
                       "EXTERNAL_POTENTIAL",error=error)
    CALL section_vals_get(ext_pot_section,explicit=explicit,error=error)
    IF (explicit) THEN
       output_unit = -1
       logger => cp_error_get_logger(error)
       ionode = logger%para_env%mepos==logger%para_env%source
       IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

       CALL initf(1)
       CALL section_vals_val_get(ext_pot_section,"FUNCTION",c_val=coupling_function,error=error)
       CALL compress(coupling_function, full=.TRUE.)
       CALL section_vals_val_get(ext_pot_section,"PARAMETERS",n_rep_val=n_rep,error=error)
       IF (n_rep /=0) THEN
          CALL section_vals_val_get(ext_pot_section,"PARAMETERS",c_vals=par,error=error)
          CALL section_vals_val_get(ext_pot_section,"VALUES",r_vals=val,error=error)
          sz_par = SIZE(par)
          CPPrecondition(SIZE(val)==sz_par,cp_failure_level,routineP,error,failure)
       ELSE
          sz_par = 0
       END IF
       CALL section_vals_val_get(ext_pot_section,"DX",r_val=dx,error=error)
       CALL section_vals_val_get(ext_pot_section,"ERROR_LIMIT",r_val=lerr,error=error)
       ALLOCATE(my_par(sz_par+3),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_val(sz_par+3),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (sz_par /= 0) THEN
          my_par(1:sz_par)=par
          my_val(1:sz_par)=val
       END IF
       my_par(1+sz_par)="X"
       my_par(2+sz_par)="Y"
       my_par(3+sz_par)="Z"
       CALL parsef(1,TRIM(coupling_function),my_par)

       ! Apply potential on all atoms, computing energy and forces
       NULLIFY(particles, subsys)
       CALL force_env_get(force_env, subsys=subsys, error=error)
       CALL cp_subsys_get(subsys, particles=particles, error=error)
       CALL force_env_get(force_env, additional_potential=energy, error=error)
       DO iatom = 1, particles%n_els
          my_val(1+sz_par)=particles%els(iatom)%r(1)
          my_val(2+sz_par)=particles%els(iatom)%r(2)
          my_val(3+sz_par)=particles%els(iatom)%r(3)
          
          energy = energy + evalf(1,my_val)
          DO j = 1, 3
             dedf = evalfd(1,sz_par+j,my_val,dx,err)
             IF ((ABS(err)>lerr).AND.(output_unit>0)) THEN
                WRITE(this_error,"(A,G12.6,A)")"(",err,")"
                WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
                CALL compress(this_error,.TRUE.)
                CALL compress(def_error,.TRUE.)
                CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
                     'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                     ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                     TRIM(def_error)//' .',error)
             END IF
             particles%els(iatom)%f(j)=particles%els(iatom)%f(j)-dedf
          END DO
       END DO
       CALL force_env_set(force_env, additional_potential=energy, error=error)
       DEALLOCATE(my_par,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_val,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL finalizef()
    END IF
    CALL timestop(handle)
  END SUBROUTINE add_external_potential

END MODULE external_potential_methods
