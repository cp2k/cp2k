!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief superstucture that hold various representations of the density and
!>      keeps track of which ones are valid
!> \par History
!>      08.2002 created [fawzi]
!>      08.2014 kpoints [JGH]
!>      11.2014 make qs_rho_type PRIVATE [Ole Schuett]
!>      11.2014 unified k-point and gamma-point code [Ole Schuett]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_rho_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE kinds,                           ONLY: dp
  USE kpoint_transitional,             ONLY: get_1d_pointer,&
                                             get_2d_pointer,&
                                             kpoint_transitional_release,&
                                             kpoint_transitional_type,&
                                             set_1d_pointer,&
                                             set_2d_pointer
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_types'

  PUBLIC :: qs_rho_p_type, qs_rho_type
  PUBLIC :: qs_rho_retain, qs_rho_release,&
            qs_rho_get, qs_rho_set, qs_rho_clear, qs_rho_create

! *****************************************************************************
!> \brief keeps the density in various representations, keeping track of
!>      which ones are valid.
!> \param most attributes are array with either lda or lsd_alpha,lsd_beta.
!> \param rho_ao the filtered rho in the localized atom basis (to have rho(r)
!>        the filtred matrix is enough, but rho(r,r') is lost).
!> \param rho_ao_kp the filtered rho in the localized atom basis (to have rho(r)
!>        the filtered matrix is enough, but rho(r,r') is lost).
!>        for kpoints, in real space index form
!> \param rho_r grids with rho in the real space
!> \param tau_r grids with the kinetic energy density in real space
!> \param rho_g grids with rho in the g space
!> \param tau_g grids with the kinetic energy density in g space
!> \param rho_g_valid , rho_r_valid, tau_r_valid, tau_g_valid: if the
!>        corresponding component is valid
!> \param ref_count the reference count, when it becomes 0 the type
!>        is deallocated.
!> \param rebuild_each how often a rebuild should be done by default
!> \param tot_rho_r the total charge in r space (valid only if rho_r is)
!> \note
!>      If pw_p_type would implement retain/release it would be nice to
!>      store also the core charge and the qs_charges in this structure...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_rho_type
     PRIVATE
     TYPE(kpoint_transitional_type)                 :: rho_ao
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER   :: rho_ao_im => Null()
     TYPE(pw_p_type), DIMENSION(:), POINTER         :: rho_g     => Null(),&
                                                       rho_r     => Null(),&
                                                       drho_g    => Null(),&
                                                       drho_r    => Null(),&
                                                       tau_g     => Null(),&
                                                       tau_r     => Null()
     ! Final rho_iter of last SCCS cycle (r-space)
     TYPE(pw_p_type), POINTER                       :: rho_r_sccs => Null()
     LOGICAL                                        :: rho_g_valid  = .FALSE., &
                                                       rho_r_valid  = .FALSE., &
                                                       drho_r_valid = .FALSE., &
                                                       drho_g_valid = .FALSE., &
                                                       tau_r_valid  = .FALSE., &
                                                       tau_g_valid  = .FALSE., &
                                                       soft_valid   = .FALSE.
     INTEGER                                        :: ref_count = -1, &
                                                       id_nr = -1, &
                                                       rebuild_each = -1
     REAL(KIND = dp), DIMENSION(:), POINTER ::         tot_rho_r => Null(),&
                                                       tot_rho_g => Null()
  END TYPE qs_rho_type

! *****************************************************************************
  TYPE qs_rho_p_type
    TYPE(qs_rho_type), POINTER         :: rho
  END TYPE qs_rho_p_type


 INTEGER, PRIVATE, SAVE :: last_rho_id_nr = 0

CONTAINS

! *****************************************************************************
!> \brief Allocates a new instance of rho.
!> \param rho ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE qs_rho_create(rho, error)
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_rho_create', &
      routineP = moduleN//':'//routineN

   IF(ASSOCIATED(rho)) STOP routineP//" rho already associated"

   ALLOCATE(rho)
   last_rho_id_nr=last_rho_id_nr+1
   rho%id_nr=last_rho_id_nr
   rho%rebuild_each=5
   rho%ref_count=1
 END SUBROUTINE qs_rho_create

! *****************************************************************************
!> \brief retains a rho_struct by increasing the reference count by one
!>      (to be called when you want to keep a shared copy)
!> \param rho_struct the structure to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_retain(rho_struct,error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_struct%ref_count=rho_struct%ref_count+1
    END IF
  END SUBROUTINE qs_rho_retain

! *****************************************************************************
!> \brief releases a rho_struct by decreasing the reference count by one
!>      and deallocating if it reaches 0 (to be called when you don't want
!>      anymore a shared copy)
!> \param rho_struct the structure to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_release(rho_struct,error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(rho_struct)) THEN
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_struct%ref_count=rho_struct%ref_count-1
       IF (rho_struct%ref_count<1) THEN
          CALL qs_rho_clear(rho_struct, error)
          DEALLOCATE(rho_struct)
       END IF
    END IF

    NULLIFY (rho_struct)

 END SUBROUTINE qs_rho_release


! *****************************************************************************
!> \brief Deallocates all components, whithout deallocating rho_struct itself.
!> \param rho_struct ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qs_rho_clear(rho_struct, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_clear', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (ASSOCIATED(rho_struct%rho_r)) THEN
       DO i=1,SIZE(rho_struct%rho_r)
          CALL pw_release(rho_struct%rho_r(i)%pw,error=error)
       END DO
       DEALLOCATE(rho_struct%rho_r)
    END IF
    IF (ASSOCIATED(rho_struct%drho_r)) THEN
       DO i=1,SIZE(rho_struct%drho_r)
          CALL pw_release(rho_struct%drho_r(i)%pw,error=error)
       END DO
       DEALLOCATE(rho_struct%drho_r)
    END IF
    IF (ASSOCIATED(rho_struct%drho_g)) THEN
       DO i=1,SIZE(rho_struct%drho_g)
          CALL pw_release(rho_struct%drho_g(i)%pw,error=error)
       END DO
       DEALLOCATE(rho_struct%drho_g)
    END IF
    IF (ASSOCIATED(rho_struct%tau_r)) THEN
       DO i=1,SIZE(rho_struct%tau_r)
          CALL pw_release(rho_struct%tau_r(i)%pw,error=error)
       END DO
       DEALLOCATE(rho_struct%tau_r)
    END IF
    IF (ASSOCIATED(rho_struct%rho_g)) THEN
       DO i=1,SIZE(rho_struct%rho_g)
          CALL pw_release(rho_struct%rho_g(i)%pw,error=error)
       END DO
       DEALLOCATE(rho_struct%rho_g)
    END IF
    IF (ASSOCIATED(rho_struct%tau_g)) THEN
       DO i=1,SIZE(rho_struct%tau_g)
          CALL pw_release(rho_struct%tau_g(i)%pw,error=error)
       END DO
       DEALLOCATE(rho_struct%tau_g)
    END IF
    IF (ASSOCIATED(rho_struct%rho_r_sccs)) THEN
       CALL pw_release(rho_struct%rho_r_sccs%pw,error=error)
       DEALLOCATE(rho_struct%rho_r_sccs)
    END IF

    CALL kpoint_transitional_release(rho_struct%rho_ao, error)

    IF (ASSOCIATED(rho_struct%rho_ao_im))&
       CALL cp_dbcsr_deallocate_matrix_set(rho_struct%rho_ao_im,error=error)
    IF (ASSOCIATED(rho_struct%tot_rho_r))&
       DEALLOCATE(rho_struct%tot_rho_r)
    IF (ASSOCIATED(rho_struct%tot_rho_g))&
       DEALLOCATE(rho_struct%tot_rho_g)

  END SUBROUTINE qs_rho_clear

! *****************************************************************************
!> \brief returns info about the density described by this object.
!>      If some representation is not available an error is issued
!> \param rho_struct ...
!> \param rho_ao ...
!> \param rho_ao_im ...
!> \param rho_ao_kp ...
!> \param rho_r ...
!> \param drho_r ...
!> \param rho_g ...
!> \param drho_g ...
!> \param tau_r ...
!> \param tau_g ...
!> \param rho_r_valid ...
!> \param drho_r_valid ...
!> \param rho_g_valid ...
!> \param drho_g_valid ...
!> \param tau_r_valid ...
!> \param tau_g_valid ...
!> \param rebuild_each ...
!> \param tot_rho_r ...
!> \param tot_rho_g ...
!> \param rho_r_sccs ...
!> \param soft_valid ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      For the other inputs the qs_rho_type attributes
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_get(rho_struct, rho_ao, rho_ao_im, rho_ao_kp, rho_r, drho_r, &
       rho_g, drho_g, tau_r, tau_g, rho_r_valid, drho_r_valid, rho_g_valid, &
       drho_g_valid, tau_r_valid, tau_g_valid, rebuild_each, tot_rho_r, tot_rho_g,&
       rho_r_sccs, soft_valid, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_ao, rho_ao_im
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rho_ao_kp
    TYPE(pw_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_r, drho_r, rho_g, drho_g, &
                                                tau_r, tau_g
    LOGICAL, INTENT(out), OPTIONAL           :: rho_r_valid, drho_r_valid, &
                                                rho_g_valid, drho_g_valid, &
                                                tau_r_valid, tau_g_valid
    INTEGER, INTENT(out), OPTIONAL           :: rebuild_each
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: tot_rho_r, tot_rho_g
    TYPE(pw_p_type), OPTIONAL, POINTER       :: rho_r_sccs
    LOGICAL, INTENT(out), OPTIONAL           :: soft_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(rho_ao))    rho_ao    => get_1d_pointer(rho_struct%rho_ao)
    IF (PRESENT(rho_ao_kp)) rho_ao_kp => get_2d_pointer(rho_struct%rho_ao)

    IF (PRESENT(rho_ao_im))    rho_ao_im    => rho_struct%rho_ao_im
    IF (PRESENT(rho_r))        rho_r        => rho_struct%rho_r
    IF (PRESENT(drho_r))       drho_r       => rho_struct%drho_r
    IF (PRESENT(rho_g))        rho_g        => rho_struct%rho_g
    IF (PRESENT(drho_g))       drho_g       => rho_struct%drho_g
    IF (PRESENT(tau_r))        tau_r        => rho_struct%tau_r
    IF (PRESENT(tau_g))        tau_g        => rho_struct%tau_g
    IF (PRESENT(rho_r_valid))  rho_r_valid  =  rho_struct%rho_r_valid
    IF (PRESENT(rho_g_valid))  rho_g_valid  =  rho_struct%rho_g_valid
    IF (PRESENT(drho_r_valid)) drho_r_valid =  rho_struct%drho_r_valid
    IF (PRESENT(drho_g_valid)) drho_g_valid =  rho_struct%drho_g_valid
    IF (PRESENT(tau_r_valid))  tau_r_valid  =  rho_struct%tau_r_valid
    IF (PRESENT(tau_g_valid))  tau_g_valid  =  rho_struct%tau_g_valid
    IF (PRESENT(soft_valid))   soft_valid   =  rho_struct%soft_valid
    IF (PRESENT(rebuild_each)) rebuild_each =  rho_struct%rebuild_each
    IF (PRESENT(tot_rho_r))    tot_rho_r    => rho_struct%tot_rho_r
    IF (PRESENT(tot_rho_g))    tot_rho_g    => rho_struct%tot_rho_g
    IF (PRESENT(rho_r_sccs))   rho_r_sccs   => rho_struct%rho_r_sccs

  END SUBROUTINE qs_rho_get

! *****************************************************************************
!> \brief ...
!> \param rho_struct ...
!> \param rho_ao ...
!> \param rho_ao_im ...
!> \param rho_ao_kp ...
!> \param rho_r ...
!> \param drho_r ...
!> \param rho_g ...
!> \param drho_g ...
!> \param tau_r ...
!> \param tau_g ...
!> \param rho_r_valid ...
!> \param drho_r_valid ...
!> \param rho_g_valid ...
!> \param drho_g_valid ...
!> \param tau_r_valid ...
!> \param tau_g_valid ...
!> \param rebuild_each ...
!> \param tot_rho_r ...
!> \param tot_rho_g ...
!> \param rho_r_sccs ...
!> \param soft_valid ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qs_rho_set(rho_struct, rho_ao, rho_ao_im, rho_ao_kp, rho_r, drho_r, &
       rho_g, drho_g, tau_r, tau_g, rho_r_valid, drho_r_valid, rho_g_valid, &
       drho_g_valid, tau_r_valid, tau_g_valid, rebuild_each, tot_rho_r, tot_rho_g,&
       rho_r_sccs, soft_valid, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_ao, rho_ao_im
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rho_ao_kp
    TYPE(pw_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_r, drho_r, rho_g, drho_g, &
                                                tau_r, tau_g
    LOGICAL, INTENT(in), OPTIONAL            :: rho_r_valid, drho_r_valid, &
                                                rho_g_valid, drho_g_valid, &
                                                tau_r_valid, tau_g_valid
    INTEGER, INTENT(in), OPTIONAL            :: rebuild_each
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: tot_rho_r, tot_rho_g
    TYPE(pw_p_type), OPTIONAL, POINTER       :: rho_r_sccs
    LOGICAL, INTENT(in), OPTIONAL            :: soft_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(rho_ao))    CALL set_1d_pointer(rho_struct%rho_ao, rho_ao)
    IF (PRESENT(rho_ao_kp)) CALL set_2d_pointer(rho_struct%rho_ao, rho_ao_kp)

    IF (PRESENT(rho_ao_im))    rho_struct%rho_ao_im    => rho_ao_im
    IF (PRESENT(rho_r))        rho_struct%rho_r        => rho_r
    IF (PRESENT(rho_g))        rho_struct%rho_g        => rho_g
    IF (PRESENT(drho_r))       rho_struct%drho_r       => drho_r
    IF (PRESENT(drho_g))       rho_struct%drho_g       => drho_g
    IF (PRESENT(tau_r))        rho_struct%tau_r        => tau_r
    IF (PRESENT(tau_g))        rho_struct%tau_g        => tau_g
    IF (PRESENT(rho_r_valid))  rho_struct%rho_r_valid  =  rho_r_valid
    IF (PRESENT(rho_g_valid))  rho_struct%rho_g_valid  =  rho_g_valid
    IF (PRESENT(drho_r_valid)) rho_struct%drho_r_valid =  drho_r_valid
    IF (PRESENT(drho_g_valid)) rho_struct%drho_g_valid =  drho_g_valid
    IF (PRESENT(tau_r_valid))  rho_struct%tau_r_valid  =  tau_r_valid
    IF (PRESENT(tau_g_valid))  rho_struct%tau_g_valid  =  tau_g_valid
    IF (PRESENT(soft_valid))   rho_struct%soft_valid   =  soft_valid
    IF (PRESENT(rebuild_each)) rho_struct%rebuild_each =  rebuild_each
    IF (PRESENT(tot_rho_r))    rho_struct%tot_rho_r    => tot_rho_r
    IF (PRESENT(tot_rho_g))    rho_struct%tot_rho_g    => tot_rho_g
    IF (PRESENT(rho_r_sccs))   rho_struct%rho_r_sccs   => rho_r_sccs

  END SUBROUTINE qs_rho_set

END MODULE qs_rho_types
