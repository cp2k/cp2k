!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Working with the  semi empirical parameter types.
!> \author JGH (14.08.2004)
! *****************************************************************************
MODULE semi_empirical_utils
  USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                             create_gto_from_sto_basis,&
                                             gto_basis_set_type,&
                                             set_sto_basis_set
  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_mndod,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_rm1
  USE kinds,                           ONLY: dp
  USE semi_empirical_mpole_methods,    ONLY: semi_empirical_mpole_p_setup
  USE semi_empirical_par_utils,        ONLY: get_se_basis,&
                                             setup_1c_2el_int
  USE semi_empirical_parameters,       ONLY: am1_default_parameter,&
                                             mndo_default_parameter,&
                                             pcharge_default_parameter,&
                                             pdg_default_parameter,&
                                             pm3_default_parameter,&
                                             pm6_default_parameter,&
                                             rm1_default_parameter
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_utils'

  PUBLIC :: init_se_param, se_param_set_default, get_se_type

CONTAINS

! *****************************************************************************
!> \brief Initialize semi_empirical type
! *****************************************************************************
  SUBROUTINE init_se_param(sep,orb_basis_set,ngauss,error)

    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    INTEGER, INTENT(IN)                      :: ngauss
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_se_param', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: symbol
    INTEGER                                  :: l, nshell, stat
    INTEGER, DIMENSION(:), POINTER           :: lq, nq
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet

    failure = .FALSE.
    IF (ASSOCIATED(sep)) THEN
       CALL allocate_sto_basis_set (sep%basis,error)
       nshell = 0
       IF (sep%natorb ==  1) nshell=1
       IF (sep%natorb ==  4) nshell=2
       IF (sep%natorb ==  9) nshell=3
       IF (sep%natorb == 16) nshell=4
       ALLOCATE (nq(0:3),lq(0:3),zet(0:3),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       ALLOCATE (symbol(0:3),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       symbol=""
       nq=0
       lq=0
       zet=0._dp
       DO l=0,nshell-1
          nq(l)=get_se_basis(sep%z,l,error)
          lq(l)=l
          zet(l)=sep%sto_exponents(l)
          IF(l==0) WRITE(symbol(0),'(I1,A1)') nq(l),"S"
          IF(l==1) WRITE(symbol(1),'(I1,A1)') nq(l),"P"
          IF(l==2) WRITE(symbol(2),'(I1,A1)') nq(l),"D"
          IF(l==3) WRITE(symbol(3),'(I1,A1)') nq(l),"F"
       END DO

       IF (nshell > 0) THEN
          sep%ngauss = ngauss
          CALL set_sto_basis_set(sep%basis,name=sep%name,nshell=nshell,symbol=symbol,&
               nq=nq,lq=lq,zet=zet)
          CALL create_gto_from_sto_basis(sep%basis,orb_basis_set,sep%ngauss,error=error)
       END IF

       DEALLOCATE (nq,lq,zet,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (symbol,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sep is not associated")
    END IF

  END SUBROUTINE init_se_param

! *****************************************************************************
!> \brief Initialize parameter for a semi_empirival type
! *****************************************************************************
  SUBROUTINE se_param_set_default(sep,z,method,error)

    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER, INTENT(IN)                      :: z, method
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_param_set_default', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(sep)) THEN
       IF (z<0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Atomic number < 0")
       END IF
       SELECT CASE (method)
       CASE (do_method_am1)
          CALL am1_default_parameter(sep,z,error)
       CASE (do_method_rm1)
          CALL rm1_default_parameter(sep,z,error)
       CASE (do_method_pm3)
          CALL pm3_default_parameter(sep,z,error)
       CASE (do_method_pm6)
          CALL pm6_default_parameter(sep,z,error)
       CASE (do_method_pdg)
          CALL pdg_default_parameter(sep,z,error)
       CASE (do_method_mndo)
          CALL mndo_default_parameter(sep,z,do_method_mndo, error)
       CASE (do_method_mndod)
          CALL mndo_default_parameter(sep,z,do_method_mndod,error)
       CASE (do_method_pchg)
          CALL pcharge_default_parameter(sep,z,error)
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Semiempirical method unknown")
       END SELECT
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sep is not associated")
    END IF

    ! Check if the element has been defined.. 
    CALL cp_assert(sep%defined,cp_fatal_level,cp_assertion_failed,routineP,&
         "Semiempirical type ("//TRIM(sep%name)//") cannot be defined for "//&
         "the requested parameterization."//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
   
    ! Fill 1 center - 2 electron integrals 
    CALL setup_1c_2el_int(sep, error)
    
    ! Fill multipolar expansion of atomic orbitals charge distributions
    CALL semi_empirical_mpole_p_setup(sep%w_mpole, sep, error)

    ! Get the value of the size of CORE integral array
    sep%core_size = 0
    IF (sep%natorb>0) sep%core_size = 1
    IF (sep%natorb>1) sep%core_size = 4
    IF (sep%dorb)     sep%core_size =10
    
    ! Get size of the all possible combinations of atomic orbitals
    sep%atm_int_size = (sep%natorb+1)*sep%natorb/2

  END SUBROUTINE se_param_set_default

! *****************************************************************************
!> \brief Gives back the unique semi_empirical METHOD type
! *****************************************************************************
  FUNCTION get_se_type(se_method) RESULT(se_type)

    INTEGER, INTENT(IN)                      :: se_method
    INTEGER                                  :: se_type

    CHARACTER(len=*), PARAMETER :: routineN = 'get_se_type', &
      routineP = moduleN//':'//routineN

    SELECT CASE (se_method)
    CASE DEFAULT
       se_type=se_method
    CASE (do_method_am1,do_method_rm1)
       se_type=do_method_am1
    END SELECT

  END FUNCTION get_se_type

END MODULE semi_empirical_utils

