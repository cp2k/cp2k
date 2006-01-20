!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_p_env_types [1.0] *
!!
!!   NAME
!!     qs_p_env_types
!!
!!   FUNCTION
!!     basis types for the calculation of the perturbation of density theory.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
! ****************************************************************************
MODULE qs_p_env_types
  USE cp_fm_types,                     ONLY: cp_fm_release,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc,&
                                             cp_fm_vect_write
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hartree_local_types,             ONLY: hartree_local_release,&
                                             hartree_local_type
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             preconditioner_type
  USE qs_kpp1_env_types,               ONLY: kpp1_release,&
                                             qs_kpp1_env_type
  USE qs_local_rho_types,              ONLY: local_rho_set_release,&
                                             local_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_type,&
                                             qs_rho_release
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qs_p_env_type, qs_p_env_p_type
  PUBLIC :: p_env_retain, p_env_release, p_env_write

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_env_types'

!!***
!****************************************************************************


!!****s* qs_p_env_types/qs_p_env_type [1.0] *
!!
!!   NAME
!!     qs_p_env_type
!!
!!   FUNCTION
!!     Represent a qs system that is perturbed.
!!     Can calculate the linear operator and the rhs of the system 
!!     of equations that needs to be solved for the perturbation.
!!
!!   NOTES
!!     for the moment no smearing of the orbitals.
!!
!!   ATTRIBUTES
!!     - orthogonal_orbitals: if the orbitals are orthogonal
!!     - id_nr: identification number (unique between p_env_types)
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - iter: number of iterations
!!     - kpp1: the kpp1 matrix (make it temporary?)
!!     - m_epsilon: minus epsilon: -<psi0d|H_rho|psi0d>
!!     - psi0d: the dual basis of psi0: psi0 (psi0^T S psi0)^-1
!!     - S_psi0: S times psi0, cached for performace reasons
!!     - Smo_inv: inverse of the mo overlap: (psi0^T S psi0)^-1
!!     - kpp1_env: environment for the calculation of kpp1
!!     - rho1: the density rho1
!!     - rho1: the soft density rho1 for gapw_xc
!!     - n_mo: cached number of mo: n_mo(i)=qs_env%c(i)%nmo
!!     - n_ao: cached number of ao: n_ao(i)=qs_env%c(i)%nao
!!
!!   SOURCE
!***************************************************************************

  TYPE qs_p_env_type

     LOGICAL :: orthogonal_orbitals
     INTEGER :: id_nr, ref_count, iter
     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: kpp1, p1
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: m_epsilon,&
          psi0d, S_psi0, Smo_inv
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
     TYPE(qs_rho_type), POINTER :: rho1
     TYPE(qs_rho_type), POINTER :: rho1_xc
     INTEGER, DIMENSION(2) :: n_mo, & ! no of molecular orbitals
                              n_ao    ! no of basis functions
     ! GAPW stuff
     TYPE(hartree_local_type),  POINTER                 :: hartree_local
     TYPE(local_rho_type), POINTER                      :: local_rho_set

     ! Linear Response Modules

     TYPE(cp_fm_p_type), DIMENSION(:), POINTER  :: PS_psi0
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER  :: d_psi1, g_psi1, g_psi1_old
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER  :: buf1, buf2, os_psi0 
     TYPE(cp_fm_type), POINTER :: sinv_test

     ! preconditioner matrix should be symmetric and positive definite
     TYPE(preconditioner_type), POINTER :: preconditioner
     LOGICAL :: new_preconditioner

     !factors
     REAL(KIND = dp)  :: delta,gnorm,gnorm_cross,gnorm_old,etotal,gradient
     !line search
     INTEGER :: ls_count
     REAL(KIND = dp)  :: ls_pos(53),ls_energy(53),ls_grad(53)
     LOGICAL :: only_energy, os_valid

  END TYPE qs_p_env_type

!!***
!****************************************************************************

!!****s* qs_p_env_types/qs_p_env_p_type [1.0] *
!!
!!   NAME
!!     qs_p_env_p_type
!!
!!   FUNCTION
!!     to have arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     p_env: the pointer to the p_env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_p_env_p_type
     TYPE(qs_p_env_type), POINTER :: p_env
  END TYPE qs_p_env_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* qs_p_env_types/p_env_retain [1.0] *
!!
!!   NAME
!!     p_env_retain
!!
!!   SYNOPSIS
!!     Subroutine p_env_retain(p_env, error)
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine p_env_retain
!!
!!   FUNCTION
!!     retains the given p_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the p_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE p_env_retain(p_env,error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
     p_env%ref_count=p_env%ref_count+1
  END IF
END SUBROUTINE p_env_retain
!***************************************************************************

!!****f* qs_p_env_types/p_env_release [1.0] *
!!
!!   NAME
!!     p_env_release
!!
!!   SYNOPSIS
!!     Subroutine p_env_release(p_env, error)
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine p_env_release
!!
!!   FUNCTION
!!     relases the given p_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE p_env_release(p_env, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, ispin, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(p_env)) THEN
       CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
       p_env%ref_count=p_env%ref_count-1
       IF (p_env%ref_count<1) THEN
          CALL kpp1_release(p_env%kpp1_env, error)
          CALL cp_fm_vect_dealloc(p_env%S_psi0, error=error)
          CALL cp_fm_vect_dealloc(p_env%m_epsilon, error=error)
          CALL cp_fm_vect_dealloc(p_env%psi0d, error=error)
          CALL cp_fm_vect_dealloc(p_env%Smo_inv, error=error)
          IF(ASSOCIATED(p_env%rho1_xc)) THEN
            NULLIFY(p_env%rho1_xc%rho_ao)
            CALL qs_rho_release(p_env%rho1_xc)
          ENDIF
          CALL qs_rho_release(p_env%rho1)
          IF (ASSOCIATED(p_env%kpp1)) CALL deallocate_matrix_set(p_env%kpp1)
          IF (ASSOCIATED(p_env%p1)) CALL deallocate_matrix_set(p_env%p1)
          IF (ASSOCIATED(p_env%local_rho_set)) THEN
             CALL local_rho_set_release(p_env%local_rho_set)
          END IF
          IF (ASSOCIATED(p_env%hartree_local)) THEN
             NULLIFY(p_env%hartree_local%vtriple) ! we must not deallocate this
             CALL hartree_local_release(p_env%hartree_local)
          END IF
          IF(ASSOCIATED(p_env%PS_psi0)) THEN
            CALL cp_fm_vect_dealloc(p_env%PS_psi0, error=error)  
          END IF
          IF(ASSOCIATED(p_env%g_psi1_old)) THEN
            CALL cp_fm_vect_dealloc(p_env%g_psi1_old, error=error)  
          END IF
          IF(ASSOCIATED(p_env%g_psi1)) THEN
            CALL cp_fm_vect_dealloc(p_env%g_psi1, error=error)  
          END IF
          IF(ASSOCIATED(p_env%d_psi1)) THEN
            CALL cp_fm_vect_dealloc(p_env%d_psi1, error=error)  
          END IF
          IF(ASSOCIATED(p_env%buf1)) THEN
            CALL cp_fm_vect_dealloc(p_env%buf1, error=error)
          END IF
          IF(ASSOCIATED(p_env%buf2)) THEN
            CALL cp_fm_vect_dealloc(p_env%buf2, error=error)
          END IF
          IF(ASSOCIATED(p_env%os_psi0)) THEN
            CALL cp_fm_vect_dealloc(p_env%os_psi0, error=error)
          END IF
          IF(ASSOCIATED(p_env%preconditioner)) THEN
            CALL destroy_preconditioner(p_env%preconditioner)
          END IF
          IF(ASSOCIATED(p_env%sinv_test)) THEN
            CALL cp_fm_release(p_env%sinv_test)
          END IF
       END IF
       DEALLOCATE(p_env,stat=stat)
       CPPostcondition((stat==0),cp_failure_level,routineP,error,failure)
    END IF
    NULLIFY(p_env)
  END SUBROUTINE p_env_release
!***************************************************************************

!!****f* qs_p_env_types/p_env_write [1.0] *
!!
!!   NAME
!!     p_env_write
!!
!!   SYNOPSIS
!!     Subroutine p_env_write(p_env, unit_nr, long_description, local,&
!!         error)
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Optional, Intent (IN):: long_description, local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine p_env_write
!!
!!   FUNCTION
!!     writes the content of the p_env to the given output unit
!!
!!   NOTES
!!     remove all iostat checks??
!!
!!   ARGUMENTS
!!     - p_env: the p_env to output
!!     - unit_nr: the unit where to output
!!     - long_description: if a long description should be written
!!       (defaults to false)
!!     - local: if the unit is a local unit or a global unit 
!!       (defaults to false, i.e. global)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE p_env_write(p_env, unit_nr, long_description, local, error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description, local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iostat
    INTEGER, SAVE                            :: uniquing_nr = 0
    LOGICAL                                  :: failure, my_local, &
                                                my_long_description, should_w
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure=.FALSE.; my_local=.FALSE.; my_long_description=.FALSE.
    iostat=0

    IF (PRESENT(local)) my_local=local
    IF (PRESENT(long_description)) my_long_description=long_description
    logger=>cp_error_get_logger(error) 
    para_env=>logger%para_env
    should_w=my_local .OR. para_env%mepos==para_env%source

    IF (ASSOCIATED(p_env)) THEN
       IF (should_w) THEN
          WRITE (unit=unit_nr,&
               fmt="(' <p_env>:{ id_nr=',i10,' ref_count=',i10,',')",&
               iostat=iostat) p_env%id_nr, p_env%ref_count
          CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
          WRITE (unit=unit_nr,&
               fmt="(' iter=',i10,', orthogonal_orbitals=',l1,',')",&
               iostat=iostat) p_env%iter,p_env%orthogonal_orbitals
          CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
          WRITE (unit=unit_nr,&
               fmt="(' n_mo=',2i10,', n_ao=',2i10,',')",&
               iostat=iostat) p_env%n_mo, p_env%n_ao
          CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
          WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " m_epsilon="
          CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
       END IF
       CALL cp_fm_vect_write(p_env%m_epsilon,unit_nr=unit_nr,&
            long_description=my_long_description, local=my_local,&
            error=error)
       IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" , psi0d="
       CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
       CALL cp_fm_vect_write(p_env%psi0d,unit_nr=unit_nr,&
            long_description=my_long_description, local=my_local,&
           error=error)
       IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" , S_psi0="
       CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
       CALL cp_fm_vect_write(p_env%S_psi0,unit_nr=unit_nr,&
            long_description=my_long_description, local=my_local,&
            error=error)
       IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" , Smo_inv="
       CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
       CALL cp_fm_vect_write(p_env%Smo_inv,unit_nr=unit_nr,&
            long_description=my_long_description, local=my_local,&
            error=error)
       IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat, advance='no')&
            " , kpp1="
       CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
       IF (ASSOCIATED(p_env%kpp1)) THEN
!FM          IF (my_long_description) THEN
!FM             IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat) " ("
!FM             CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
!FM             uniquing_nr=uniquing_nr+1
!FM             DO i=1,SIZE(p_env%kpp1)
!FM                CALL cp_sm_output(logger,&
!FM                     outputName="p_env"//cp_to_string(p_env%id_nr)//'kpp1-'//&
!FM                     cp_to_string(i),&
!FM                     fromWhere=routineP,iter=uniquing_nr,&
!FM                     matrix=p_env%kpp1(i)%matrix,&
!FM                     para_env=para_env, comment='kpp1 matrix',error=error)
!FM                IF (should_w.and.i/=SIZE(p_env%kpp1)) THEN
!FM                   WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" ,"
!FM                   CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
!FM                END IF
!FM             END DO
!FM             IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat) " ),"
!FM             CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
!FM          ELSE
             IF (should_w) THEN
                WRITE(unit=unit_nr, fmt="(a)") '*associated*,'
             END IF
!FM          END IF
       ELSE
          WRITE(unit=unit_nr, fmt="(a)") '*null*,'
       END IF
       IF (should_w) THEN
          WRITE(unit=unit_nr,fmt="(a)",iostat=iostat, advance='no')" kpp1_env="
          CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
          IF (ASSOCIATED(p_env%kpp1_env)) THEN
             WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" *associated*,"
          ELSE
             WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" *null*,"
          END IF
!FM       call kpp1_env_write(p_env%kpp1_env,unit_nr=unit_nr,&
!FM            long_description=my_long_description, local=my_local,&
!FM            error=error)
       END IF
       IF (should_w) WRITE(unit=unit_nr,fmt="(a)",iostat=iostat)" }"
       CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
    ELSE
       IF (should_w) THEN
          WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " <p_env>:*null*"
          CPInvariantNoFail(iostat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE p_env_write
!***************************************************************************

END MODULE qs_p_env_types
