!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_rho0_types
!!
!!   NAME
!!     qs_rho0_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE qs_local_rho_types

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE erf_fn,                          ONLY: erf
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE mathconstants,                   ONLY: fourpi,pi
  USE memory_utilities,                ONLY: reallocate
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_rho_atom_types,               ONLY: deallocate_rho_atom_set, &
                                             rho_atom_type
  USE qs_rho0_types,                   ONLY: deallocate_rho0_atom,&
                                             deallocate_rho0_mpole,&
                                             rho0_atom_type,&
                                             rho0_mpole_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_local_rho_types"

! *** Define rhoz and local_rho types ***

  TYPE  rhoz_type
    REAL(dp)                             ::  one_atom
    REAL(dp), DIMENSION(:), POINTER      ::  r_coef
    REAL(dp), DIMENSION(:), POINTER      ::  dr_coef
    REAL(dp), DIMENSION(:), POINTER      ::  vr_coef
  END TYPE rhoz_type

  TYPE local_rho_type
    TYPE(rho_atom_type), DIMENSION(:), POINTER            :: rho_atom_set
    TYPE(rho0_mpole_type), POINTER                        :: rho0_mpole
    TYPE(rho0_mpole_type), POINTER                        :: nddo_mpole
    TYPE(rho0_atom_type), DIMENSION(:), POINTER           :: rho0_atom_set
    TYPE(rhoz_type),  DIMENSION(:), POINTER               :: rhoz_set
    REAL(dp)                                              :: rhoz_tot
  END TYPE local_rho_type

! Public Types

  PUBLIC ::  local_rho_type, rhoz_type

! Public Subroutine
 
  PUBLIC :: allocate_rhoz, calculate_rhoz, deallocate_rhoz,& 
            get_local_rho, local_rho_set_create, &
            local_rho_set_release, set_local_rho

!******************************************************************************
 
  CONTAINS
 
!******************************************************************************
 
  SUBROUTINE allocate_rhoz(rhoz_set,nkind)

    TYPE(rhoz_type), DIMENSION(:), &
                                   POINTER :: rhoz_set
    INTEGER                                :: nkind
    
    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_rhoz"

    INTEGER                                :: ikind, istat

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(rhoz_set)) THEN
      CALL deallocate_rhoz(rhoz_set)
    END IF 

    ALLOCATE (rhoz_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "rhoz_set",nkind*int_size)

    DO ikind = 1,nkind 
      NULLIFY(rhoz_set(ikind)%r_coef) 
      NULLIFY(rhoz_set(ikind)%dr_coef)
      NULLIFY(rhoz_set(ikind)%vr_coef)
    ENDDO

  END SUBROUTINE allocate_rhoz

! ******************************************************************************

  SUBROUTINE calculate_rhoz(rhoz,grid_atom,alpha,zeff,natom,rhoz_tot,harmonics)
 
    TYPE(rhoz_type)             :: rhoz
    TYPE(grid_atom_type)        :: grid_atom
    TYPE(harmonics_atom_type)   :: harmonics
    REAL(dp), INTENT(IN)        :: alpha
    REAL(dp), INTENT(INOUT)     :: rhoz_tot
    INTEGER                     :: natom,zeff
 
    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_rhoz"

    INTEGER                     :: ir,nr,na,ia
    REAL(dp)                    :: c1, c2,c3, prefactor1, prefactor2, prefactor3
    REAL(dp)                    :: sum,t
!   ---------------------------------------------------------------------------

    nr = grid_atom%nr
    na = grid_atom%ng_sphere
    CALL reallocate(rhoz%r_coef,1,nr) 
    CALL reallocate(rhoz%dr_coef,1,nr) 
    CALL reallocate(rhoz%vr_coef,1,nr) 

    c1 = alpha/pi
    c2 = c1*c1*c1*fourpi
    c3 = SQRT(alpha)
    prefactor1 = zeff*SQRT(c2)
    prefactor2 = -2.0_dp*alpha
    prefactor3 = -zeff*SQRT(fourpi)

    sum = 0.0_dp
    DO ir = 1,nr
      c1 = -alpha*grid_atom%rad2(ir)
      rhoz%r_coef(ir) = -exp(c1)*prefactor1
      IF(ABS(rhoz%r_coef(ir)) < 1.0E-30_dp) THEN
        rhoz%r_coef(ir) = 0.0_dp
        rhoz%dr_coef(ir) = 0.0_dp
      ELSE
        rhoz%dr_coef(ir) = prefactor2*rhoz%r_coef(ir)
      END IF
      rhoz%vr_coef(ir) = prefactor3*erf(grid_atom%rad(ir)*c3)/grid_atom%rad(ir)
      sum = sum + rhoz%r_coef(ir)*grid_atom%wr(ir)
    END DO
    rhoz%one_atom=sum*harmonics%slm_int(1)
    rhoz_tot = rhoz_tot + natom*rhoz%one_atom

! WRITE(10,*) zeff, alpha
! DO ir =1,nr
!    WRITE(10,'(I5,3(1E30.16))') ir,grid_atom%rad(ir),rhoz%vr_coef(ir),rhoz%r_coef(ir)
! END DO

! STOP 'erf'
  END SUBROUTINE calculate_rhoz

! ******************************************************************************

  SUBROUTINE deallocate_rhoz(rhoz_set)

    TYPE(rhoz_type), DIMENSION(:), POINTER :: rhoz_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_rhoz"

    INTEGER                     :: ikind, istat, nkind

!   ---------------------------------------------------------------------------

    nkind = SIZE(rhoz_set)

    DO ikind = 1,nkind
      DEALLOCATE(rhoz_set(ikind)%r_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "rhoz_set%r_coef")
      DEALLOCATE(rhoz_set(ikind)%dr_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "rhoz_set%dr_coef")
      DEALLOCATE(rhoz_set(ikind)%vr_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "rhoz_set%vr_coef")
    END DO

    DEALLOCATE(rhoz_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                     "rhoz_set")
          
  END SUBROUTINE  deallocate_rhoz

! ******************************************************************************
  SUBROUTINE  get_local_rho(local_rho_set,rho_atom_set,rho0_atom_set,&
                            rho0_mpole,nddo_mpole,rhoz_set)

    TYPE(local_rho_type), INTENT(IN), TARGET     :: local_rho_set
    TYPE(rho_atom_type), DIMENSION(:),  &
                             OPTIONAL, POINTER   :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:), &
                              OPTIONAL, POINTER  :: rho0_atom_set
    TYPE(rho0_mpole_type),     OPTIONAL, POINTER :: rho0_mpole
    TYPE(rho0_mpole_type),     OPTIONAL, POINTER :: nddo_mpole
    TYPE(rhoz_type), DIMENSION(:), &
                         OPTIONAL, POINTER       :: rhoz_set

    CHARACTER(len=*), PARAMETER :: routineN = 'get_local_rho'

!   ---------------------------------------------------------------------------

    IF (PRESENT(rho_atom_set)) rho_atom_set => local_rho_set%rho_atom_set
    IF (PRESENT(rho0_atom_set)) rho0_atom_set => local_rho_set%rho0_atom_set
    IF (PRESENT(rho0_mpole))  rho0_mpole =>  local_rho_set%rho0_mpole
    IF (PRESENT(nddo_mpole))  nddo_mpole =>  local_rho_set%nddo_mpole
    IF (PRESENT(rhoz_set))  rhoz_set => local_rho_set%rhoz_set

  END SUBROUTINE get_local_rho

! ******************************************************************************

  SUBROUTINE local_rho_set_create(local_rho_set, error)

   TYPE(local_rho_type), POINTER     :: local_rho_set
       TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "local_rho_set_create" 

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------
    failure=.FALSE.

    ALLOCATE(local_rho_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineN,error,failure)
  
    NULLIFY (local_rho_set%rho_atom_set)
    NULLIFY (local_rho_set%rho0_atom_set)
    NULLIFY (local_rho_set%rho0_mpole)
    NULLIFY (local_rho_set%nddo_mpole)
    NULLIFY (local_rho_set%rhoz_set)

  END SUBROUTINE local_rho_set_create

! ******************************************************************************
  SUBROUTINE local_rho_set_release(local_rho_set,error)

   TYPE(local_rho_type), POINTER         :: local_rho_set
   TYPE(cp_error_type), INTENT(inout), &
                      OPTIONAL           :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'local_rho_set_release', &
      routineP = module_name//':'//routineN

    INTEGER                                  :: stat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(local_rho_set)) THEN
      IF (ASSOCIATED(local_rho_set%rho_atom_set)) THEN
         CALL deallocate_rho_atom_set(local_rho_set%rho_atom_set)
      END IF

      IF (ASSOCIATED(local_rho_set%rho0_atom_set)) THEN
        CALL deallocate_rho0_atom(local_rho_set%rho0_atom_set)
      END IF

      IF (ASSOCIATED(local_rho_set%rho0_mpole)) THEN
        CALL deallocate_rho0_mpole(local_rho_set%rho0_mpole)
      END IF

      IF (ASSOCIATED(local_rho_set%nddo_mpole)) THEN
        CALL deallocate_rho0_mpole(local_rho_set%nddo_mpole)
      END IF

      IF(ASSOCIATED(local_rho_set%rhoz_set)) THEN
        CALL deallocate_rhoz(local_rho_set%rhoz_set)
      ENDIF

      DEALLOCATE(local_rho_set,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF 

  END SUBROUTINE local_rho_set_release

! ******************************************************************************

  SUBROUTINE set_local_rho(local_rho_set,rho_atom_set,rho0_atom_set,&
                            rho0_mpole,nddo_mpole,rhoz_set)

    TYPE(local_rho_type), INTENT(INout), TARGET     :: local_rho_set
    TYPE(rho_atom_type), DIMENSION(:),  &
                                OPTIONAL, POINTER   :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:), &
                                 OPTIONAL, POINTER  :: rho0_atom_set
    TYPE(rho0_mpole_type),        OPTIONAL, POINTER :: rho0_mpole
    TYPE(rho0_mpole_type),        OPTIONAL, POINTER :: nddo_mpole
    TYPE(rhoz_type),  DIMENSION(:), &
                            OPTIONAL, POINTER       :: rhoz_set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_local_rho'

!   ---------------------------------------------------------------------------

    IF (PRESENT(rho_atom_set)) THEN
       IF(ASSOCIATED(local_rho_set%rho_atom_set)) THEN
         CALL deallocate_rho_atom_set(local_rho_set%rho_atom_set)
       ENDIF
       local_rho_set%rho_atom_set => rho_atom_set
    END IF   

    IF (PRESENT(rho0_atom_set)) THEN
       IF(ASSOCIATED(local_rho_set%rho0_atom_set)) THEN
         CALL deallocate_rho0_atom(local_rho_set%rho0_atom_set)
       ENDIF
       local_rho_set%rho0_atom_set => rho0_atom_set
    END IF 

    IF (PRESENT(rho0_mpole)) THEN
       IF(ASSOCIATED(local_rho_set%rho0_mpole)) THEN
         CALL deallocate_rho0_mpole(local_rho_set%rho0_mpole)
       ENDIF
       local_rho_set%rho0_mpole => rho0_mpole
    END IF 

    IF (PRESENT(nddo_mpole)) THEN
       IF(ASSOCIATED(local_rho_set%nddo_mpole)) THEN
         CALL deallocate_rho0_mpole(local_rho_set%nddo_mpole)
       ENDIF
       local_rho_set%nddo_mpole => nddo_mpole
    END IF 

    IF (PRESENT(rhoz_set)) THEN
       IF(ASSOCIATED(local_rho_set%rhoz_set)) THEN
         CALL deallocate_rhoz(local_rho_set%rhoz_set)
       ENDIF
       local_rho_set%rhoz_set => rhoz_set
    END IF

  END SUBROUTINE set_local_rho

! ******************************************************************************

END MODULE qs_local_rho_types

!******************************************************************************
