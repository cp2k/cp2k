!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief The Harris force , which contains all the informations needed to perform
!>      and analyse the Harris functional force correction.
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE harris_force_types

  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'harris_force_types'
  INTEGER, PRIVATE, SAVE :: last_harris_id=0

  ! *** Public data types ***
  PUBLIC :: harris_force_type

  ! *** Public subroutines ***
  PUBLIC :: harris_force_release,&
            harris_force_clear,&
            harris_force_create

! *****************************************************************************
!> \brief The Harris force type
!> \param id_nr unique identifier
!> \param ref_count reference counter (see doc/ReferenceCounting.html)
!> \param f_trace The forces of the sum of eigenvalue term
!> \param f_hartree The forces due to the hartree repulsion of the electrons
!> \param f_integral_vxc The forces due to the integral over the vxc potential
!> \param f_delta_integral_vxc The differnece between f_integral_vxc & f_xc
!> \param f_xc The forces due to exchange & correlation
!> \param f_ovrl The forces due to the core_overlap
!> \param f_self The forces dur to the self-interaction
!> \param f_rho_core The forces due to the hartree repulsion of the cores
!> \param EII The forces due to the ion-ion repulsion
!> \param f_harris The total forces due to the harris functional
!> \param f_cross_integrate_v The forces of the integral of the cross terms
!> \param f_cross_integrate_v_core The forces of the integral of the cross terms
!> \param f_total The sum of the total forces of all involved atoms
!>          In principle this sum should be zero, due to translation invariance
!>          It can therefore be used as a criterion for the quality of the forces
!> \param f_kinetic The forces due to the kinetic energy
!> \param f_gth_pp The forces of the GTH pseudopotential
!> \param f_overlap The forces due to the overlap
!> \param f_nsc The non_self_consistent forces
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  TYPE harris_force_type
    INTEGER                                          :: id_nr, &
                                                        ref_count
    REAL(KIND=dp), DIMENSION(:,:), POINTER           :: f_trace, &
                                                        f_V, &
                                                        f_hartree, &
                                                        f_integral_vxc, &
                                                        f_delta_integral_vxc, &
                                                        f_xc, &
                                                        f_ovrl, &
                                                        f_self, &
                                                        f_rho_core, &
                                                        f_EII, &
                                                        f_harris, &
                                                        f_cross_integrate_v, &
                                                        f_cross_integrate_v_core, &
                                                        f_kinetic, &
                                                        f_gth_pp, &
                                                        f_overlap, &
                                                        f_nsc
    REAL(KIND=dp), DIMENSION(:), POINTER             :: f_total

  END TYPE harris_force_type

! *****************************************************************************
!> \brief Type which allows the use of arrays of pointers
!> \param harris_force_ptr Pointer to the actual harris force type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  TYPE harris_force_p_type
     TYPE(harris_force_type), POINTER :: harris_force_ptr
  END TYPE harris_force_p_type

CONTAINS

! *****************************************************************************
!> \brief Retains a harris force type (see doc/ReferenceCounting.html)
!> \param harris_force The harris force type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_force_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_retain(harris_force, error)

    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
      harris_force%ref_count = harris_force%ref_count+1
    END IF
  END SUBROUTINE harris_force_retain

! *****************************************************************************
!> \brief Releases the given harris force type (see doc/ReferenceCounting.html)
!> \param harris_force The harris force type to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_force_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_release(harris_force, error)

    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    IF (ASSOCIATED(harris_force)) THEN
      CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
      harris_force%ref_count = harris_force%ref_count-1
      IF (harris_force%ref_count<1) THEN
        IF (ASSOCIATED(harris_force%f_trace)) THEN
          DEALLOCATE(harris_force%f_trace, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_V)) THEN
          DEALLOCATE(harris_force%f_V, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_hartree)) THEN
          DEALLOCATE(harris_force%f_hartree, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_integral_vxc)) THEN
          DEALLOCATE(harris_force%f_integral_vxc, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_delta_integral_vxc)) THEN
          DEALLOCATE(harris_force%f_delta_integral_vxc, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_xc)) THEN
          DEALLOCATE(harris_force%f_xc, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_ovrl)) THEN
          DEALLOCATE(harris_force%f_ovrl, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_self)) THEN
          DEALLOCATE(harris_force%f_self, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_rho_core)) THEN
          DEALLOCATE(harris_force%f_rho_core, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_EII)) THEN
          DEALLOCATE(harris_force%f_EII, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_harris)) THEN
          DEALLOCATE(harris_force%f_harris, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_cross_integrate_v)) THEN
          DEALLOCATE(harris_force%f_cross_integrate_v, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_cross_integrate_v_core)) THEN
          DEALLOCATE(harris_force%f_cross_integrate_v_core, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_total)) THEN
          DEALLOCATE(harris_force%f_total, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_kinetic)) THEN
          DEALLOCATE(harris_force%f_kinetic, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_gth_pp)) THEN
          DEALLOCATE(harris_force%f_gth_pp, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_overlap)) THEN
          DEALLOCATE(harris_force%f_overlap, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(harris_force%f_nsc)) THEN
          DEALLOCATE(harris_force%f_nsc, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        DEALLOCATE(harris_force, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(harris_force)

  END SUBROUTINE harris_force_release

! *****************************************************************************
!> \brief Returns various attributes of the harris force type
!> \param harris_force The enquired harris force type
!> \param id_nr unique identifier
!> \param f_trace The forces of the sum of eigenvalue term
!> \param f_V ...
!> \param f_hartree The forces due to the hartree repulsion of the electrons
!> \param f_integral_vxc The forces due to the integral over the vxc potential
!> \param f_delta_integral_vxc The differnece between f_integral_vxc & f_xc
!> \param f_xc The forces due to exchange & correlation
!> \param f_ovrl The forces due to the overalp
!> \param f_self The forces dur to the self-interaction
!> \param f_rho_core The forces due to the hartree repulsion of the cores
!> \param f_EII The forces due to the ion-ion repulsion
!> \param f_harris The total forces due to the harris functional
!> \param f_cross_integrate_v The forces of the integral of the cross terms
!> \param f_cross_integrate_v_core The forces of the integral of the cross terms
!> \param f_kinetic The forces due to the kinetic energy
!> \param f_gth_pp The forces of the GTH pseudopotential
!> \param f_overlap The forces due to the overlap
!> \param f_nsc The non-self-consistent forces
!>
!>      For possible missing arguments see the attributes of harris_force_type
!> \param f_total The sum of the total forces of all involved atoms
!>          In principle this sum should be zero, due to translation invariance
!>          It can therefore be used as a criterion for the quality of the forces
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_get(harris_force, id_nr, f_trace, f_V, f_hartree, &
                              f_integral_vxc, f_delta_integral_vxc, f_xc, &
                              f_ovrl, f_self, f_rho_core, f_EII, f_harris, &
                              f_cross_integrate_v, f_cross_integrate_v_core, &
                              f_kinetic, f_gth_pp, f_overlap, f_nsc, f_total, &
                              error)

    TYPE(harris_force_type), POINTER         :: harris_force
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER :: f_trace, f_V, &
      f_hartree, f_integral_vxc, f_delta_integral_vxc, f_xc, f_ovrl, f_self, &
      f_rho_core, f_EII, f_harris, f_cross_integrate_v, &
      f_cross_integrate_v_core, f_kinetic, f_gth_pp, f_overlap, f_nsc
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: f_total
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = harris_force%id_nr
      IF (PRESENT(f_trace)) f_trace = harris_force%f_trace
      IF (PRESENT(f_V)) f_V = harris_force%f_V
      IF (PRESENT(f_hartree)) f_hartree = harris_force%f_hartree
      IF (PRESENT(f_integral_vxc)) f_integral_vxc = harris_force%f_integral_vxc
      IF (PRESENT(f_delta_integral_vxc)) f_delta_integral_vxc = &
           harris_force%f_delta_integral_vxc
      IF (PRESENT(f_xc)) f_xc = harris_force%f_xc
      IF (PRESENT(f_ovrl)) f_ovrl = harris_force%f_ovrl
      IF (PRESENT(f_self)) f_self = harris_force%f_self
      IF (PRESENT(f_rho_core)) f_rho_core = harris_force%f_rho_core
      IF (PRESENT(f_EII)) f_EII = harris_force%f_EII
      IF (PRESENT(f_harris)) f_harris = harris_force%f_harris
      IF (PRESENT(f_total)) f_total = harris_force%f_total
      IF (PRESENT(f_cross_integrate_v)) f_cross_integrate_v = &
           harris_force%f_cross_integrate_v
      IF (PRESENT(f_cross_integrate_v_core)) f_cross_integrate_v_core = &
           harris_force%f_cross_integrate_v_core
      IF (PRESENT(f_kinetic)) f_kinetic = harris_force%f_kinetic
      IF (PRESENT(f_gth_pp)) f_gth_pp = harris_force%f_gth_pp
      IF (PRESENT(f_overlap)) f_overlap = harris_force%f_overlap
      IF (PRESENT(f_nsc)) f_nsc = harris_force%f_nsc
    END IF
  END SUBROUTINE harris_force_get

! *****************************************************************************
!> \brief Sets various attributes of the harris force type
!> \param harris_force The enquired harris force type
!> \param id_nr unique identifier
!> \param f_trace The forces of the sum of eigenvalue term
!> \param f_V ...
!> \param f_hartree The forces due to the hartree repulsion of the electrons
!> \param f_integral_vxc The forces due to the integral over the vxc potential
!> \param f_delta_integral_vxc The differnece between f_integral_vxc & f_xc
!> \param f_xc The forces due to exchange & correlation
!> \param f_ovrl The forces due to the overalp
!> \param f_self The forces dur to the self-interaction
!> \param f_rho_core The forces due to the hartree repulsion of the cores
!> \param f_EII The forces due to the ion-ion repulsion
!> \param f_harris The total forces due to the harris functional
!> \param f_cross_integrate_v The forces of the integral of the cross terms
!> \param f_cross_integrate_v_core The forces of the integral of the cross terms
!> \param f_kinetic The forces due to the kinetic energy
!> \param f_gth_pp The forces of the GTH pseudopotential
!> \param f_overlap The forces due to the overlap
!> \param f_nsc The non-self-consistent forces
!>
!>      For possible missing arguments see the attributes of harris_force_type
!> \param f_total The sum of the total forces of all involved atoms
!>          In principle this sum should be zero, due to translation invariance
!>          It can therefore be used as a criterion for the quality of the forces
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_set(harris_force, id_nr, f_trace, f_V, f_hartree, &
                              f_integral_vxc, f_delta_integral_vxc, f_xc, &
                              f_ovrl, f_self, f_rho_core, f_EII, f_harris, &
                              f_cross_integrate_v, f_cross_integrate_v_core, &
                              f_kinetic, f_gth_pp, f_overlap, f_nsc, f_total, &
                              error)

    TYPE(harris_force_type), POINTER         :: harris_force
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER :: f_trace, f_V, &
      f_hartree, f_integral_vxc, f_delta_integral_vxc, f_xc, f_ovrl, f_self, &
      f_rho_core, f_EII, f_harris, f_cross_integrate_v, &
      f_cross_integrate_v_core, f_kinetic, f_gth_pp, f_overlap, f_nsc
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: f_total
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) harris_force%id_nr = id_nr
      IF (PRESENT(f_trace)) harris_force%f_trace = f_trace
      IF (PRESENT(f_V)) harris_force%f_V = f_V
      IF (PRESENT(f_hartree)) harris_force%f_hartree = f_hartree
      IF (PRESENT(f_integral_vxc)) harris_force%f_integral_vxc = f_integral_vxc
      IF (PRESENT(f_delta_integral_vxc)) harris_force%f_delta_integral_vxc &
           = f_delta_integral_vxc
      IF (PRESENT(f_xc)) harris_force%f_xc = f_xc
      IF (PRESENT(f_ovrl)) harris_force%f_ovrl = f_ovrl
      IF (PRESENT(f_self)) harris_force%f_self = f_self
      IF (PRESENT(f_rho_core)) harris_force%f_rho_core = f_rho_core
      IF (PRESENT(f_EII)) harris_force%f_EII = f_EII
      IF (PRESENT(f_harris)) harris_force%f_harris = f_harris
      IF (PRESENT(f_total)) harris_force%f_total = f_total
      IF (PRESENT(f_cross_integrate_v)) harris_force%f_cross_integrate_v &
           = f_cross_integrate_v
      IF (PRESENT(f_cross_integrate_v_core)) harris_force%f_cross_integrate_v_core &
           = f_cross_integrate_v_core
      IF (PRESENT(f_kinetic)) harris_force%f_kinetic = f_kinetic
      IF (PRESENT(f_gth_pp)) harris_force%f_gth_pp = f_gth_pp
      IF (PRESENT(f_overlap)) harris_force%f_overlap = f_overlap
      IF (PRESENT(f_nsc)) harris_force%f_nsc = f_nsc
    END IF
  END SUBROUTINE harris_force_set

! *****************************************************************************
!> \brief Reinitializes the harris force type
!> \param harris_force The harris force type to be cleared
!> \param natom The number of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_force_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_clear(harris_force, natom, error)

    TYPE(harris_force_type), POINTER         :: harris_force
    INTEGER                                  :: natom
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_clear', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      DO i = 1, natom
        DO j = 1, 3
          harris_force%f_trace(i,j)                  = 0.0_dp
          harris_force%f_V(i,j)                      = 0.0_dp
          harris_force%f_hartree(i,j)                = 0.0_dp
          harris_force%f_integral_vxc(i,j)           = 0.0_dp
          harris_force%f_delta_integral_vxc(i,j)     = 0.0_dp
          harris_force%f_xc(i,j)                     = 0.0_dp
          harris_force%f_ovrl(i,j)                   = 0.0_dp
          harris_force%f_self(i,j)                   = 0.0_dp
          harris_force%f_rho_core(i,j)               = 0.0_dp
          harris_force%f_EII(i,j)                    = 0.0_dp
          harris_force%f_harris(i,j)                 = 0.0_dp
          harris_force%f_cross_integrate_v(i,j)      = 0.0_dp
          harris_force%f_cross_integrate_v_core(i,j) = 0.0_dp
          harris_force%f_kinetic(i,j)                = 0.0_dp
          harris_force%f_gth_pp(i,j)                 = 0.0_dp
          harris_force%f_overlap(i,j)                = 0.0_dp
          harris_force%f_nsc(i,j)                    = 0.0_dp
        END DO
      END DO
      harris_force%f_total(:) = 0.0_dp
    ELSE
      CALL harris_force_create(harris_force, natom, error=error)
    END IF
  END SUBROUTINE harris_force_clear

! *****************************************************************************
!> \brief Creates a harris force type instance
!> \param harris_force The harris force type to be created
!> \param natom The number of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_force_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_create(harris_force, natom, error)

    TYPE(harris_force_type), POINTER         :: harris_force
    INTEGER                                  :: natom
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    ALLOCATE(harris_force, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_trace(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_V(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_hartree(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_integral_vxc(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_delta_integral_vxc(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_xc(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_ovrl(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_self(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_rho_core(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_EII(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_harris(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_cross_integrate_v(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_cross_integrate_v_core(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_total(3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_kinetic(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_gth_pp(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_overlap(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(harris_force%f_nsc(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      harris_force%ref_count = 1
      last_harris_id=last_harris_id+1
      harris_force%id_nr = last_harris_id
      CALL harris_force_clear(harris_force, natom, error=error)
    END IF
  END SUBROUTINE harris_force_create

END MODULE harris_force_types
