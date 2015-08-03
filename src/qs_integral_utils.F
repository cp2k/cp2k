!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Some utility functions for the calculation of integrals
!> \par History
!>      JGH: initial version 
!> \author JGH (10.07.2014)
! *****************************************************************************
MODULE qs_integral_utils

  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integral_utils'

! *** Interfaces ***

  INTERFACE get_memory_usage
     MODULE PROCEDURE get_memory_usage_a, get_memory_usage_ab, &
                      get_memory_usage_abc, get_memory_usage_abcd
  END INTERFACE

! *** Public subroutines ***

  PUBLIC :: get_memory_usage, basis_set_list_setup

CONTAINS

! *****************************************************************************
!> \brief Return the maximum memory usage in integral calculations
!> \param qs_kind_set The info for all atomic kinds
!> \param basis_type_a  Type of basis
!> \param error ...
!> \retval ldmem          Result
! *****************************************************************************
  FUNCTION get_memory_usage_a(qs_kind_set,basis_type_a,error) RESULT(ldmem)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    CHARACTER(LEN=*), INTENT(IN)             :: basis_type_a
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: ldmem

    CHARACTER(len=*), PARAMETER :: routineN = 'get_memory_usage_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: maxc, maxl, maxs

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=maxc, maxlgto=maxl, maxsgf=maxs,&
                             basis_type=basis_type_a,error=error)
    ldmem = MAX(maxc,maxs)

    CALL init_orbital_pointers(maxl+2)

  END FUNCTION get_memory_usage_a

! *****************************************************************************
!> \brief Return the maximum memory usage in integral calculations
!> \param qs_kind_set The info for all atomic kinds
!> \param basis_type_a  Type of basis
!> \param basis_type_b  Type of basis
!> \param error ...
!> \retval ldmem          Result
! *****************************************************************************
  FUNCTION get_memory_usage_ab(qs_kind_set,basis_type_a,basis_type_b,error) RESULT(ldmem)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    CHARACTER(LEN=*), INTENT(IN)             :: basis_type_a, basis_type_b
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: ldmem

    CHARACTER(len=*), PARAMETER :: routineN = 'get_memory_usage_ab', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lda, ldb

    lda = get_memory_usage_a(qs_kind_set,basis_type_a,error)
    ldb = get_memory_usage_a(qs_kind_set,basis_type_b,error)
    ldmem = MAX(lda,ldb)

  END FUNCTION get_memory_usage_ab

! *****************************************************************************
!> \brief Return the maximum memory usage in integral calculations
!> \param qs_kind_set The info for all atomic kinds
!> \param basis_type_a  Type of basis
!> \param basis_type_b  Type of basis
!> \param basis_type_c  Type of basis
!> \param error ...
!> \retval ldmem          Result
! *****************************************************************************
  FUNCTION get_memory_usage_abc(qs_kind_set,basis_type_a,&
                                basis_type_b,basis_type_c,error) RESULT(ldmem)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    CHARACTER(LEN=*), INTENT(IN)             :: basis_type_a, basis_type_b, &
                                                basis_type_c
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: ldmem

    CHARACTER(len=*), PARAMETER :: routineN = 'get_memory_usage_abc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lda, ldb, ldc

    lda = get_memory_usage_a(qs_kind_set,basis_type_a,error)
    ldb = get_memory_usage_a(qs_kind_set,basis_type_b,error)
    ldc = get_memory_usage_a(qs_kind_set,basis_type_c,error)
    ldmem = MAX(lda,ldb,ldc)

  END FUNCTION get_memory_usage_abc

! *****************************************************************************
!> \brief Return the maximum memory usage in integral calculations
!> \param qs_kind_set The info for all atomic kinds
!> \param basis_type_a  Type of basis
!> \param basis_type_b  Type of basis
!> \param basis_type_c  Type of basis
!> \param basis_type_d  Type of basis
!> \param error ...
!> \retval ldmem          Result
! *****************************************************************************
  FUNCTION get_memory_usage_abcd(qs_kind_set,basis_type_a,&
           basis_type_b,basis_type_c,basis_type_d,error) RESULT(ldmem)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    CHARACTER(LEN=*), INTENT(IN)             :: basis_type_a, basis_type_b, &
                                                basis_type_c, basis_type_d
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: ldmem

    CHARACTER(len=*), PARAMETER :: routineN = 'get_memory_usage_abcd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lda, ldb, ldc, ldd

    lda = get_memory_usage_a(qs_kind_set,basis_type_a,error)
    ldb = get_memory_usage_a(qs_kind_set,basis_type_b,error)
    ldc = get_memory_usage_a(qs_kind_set,basis_type_c,error)
    ldd = get_memory_usage_a(qs_kind_set,basis_type_d,error)
    ldmem = MAX(lda,ldb,ldc,ldd)

  END FUNCTION get_memory_usage_abcd

! *****************************************************************************

! *****************************************************************************
!> \brief RSet up an easy accessible list of the basis sets for all kinds
!> \param basis_set_list    The basis set list
!> \param basis_type ...
!> \param qs_kind_set   Kind information, the basis is used
!> \param error ...
! *****************************************************************************
  SUBROUTINE basis_set_list_setup(basis_set_list,basis_type,qs_kind_set,error)

    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    CHARACTER(len=*), INTENT(IN)             :: basis_type
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_set_list_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

! set up basis sets

    DO ikind=1,SIZE(qs_kind_set)
       qs_kind => qs_kind_set(ikind)
       CALL get_qs_kind(qs_kind=qs_kind,basis_set=basis_set,&
                        basis_type=basis_type,error=error)
       NULLIFY(basis_set_list(ikind)%gto_basis_set)
       IF (ASSOCIATED(basis_set)) basis_set_list(ikind)%gto_basis_set => basis_set
    END DO

  END SUBROUTINE basis_set_list_setup

! *****************************************************************************

END MODULE qs_integral_utils

