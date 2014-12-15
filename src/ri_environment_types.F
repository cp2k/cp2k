!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the types and subroutines for dealing with the ri_env
!>        ri : resolution of the identity
!> \author JGH July 2009
! *****************************************************************************
MODULE ri_environment_types
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_release,&
                                             semi_empirical_si_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE ri_kind_vector_type
     INTEGER                                          :: natom
     INTEGER                                          :: nbasis
     REAL(dp), DIMENSION(:,:), POINTER                :: vmat
  END TYPE ri_kind_vector_type

  TYPE ri_vector_type
     INTEGER                                          :: total_length
     TYPE(ri_kind_vector_type), DIMENSION(:), POINTER :: vector
  END TYPE ri_vector_type

  TYPE ri_environment_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE(ri_vector_type), POINTER                    :: coeff
     TYPE(ri_vector_type), POINTER                    :: rhs
     REAL(KIND=dp)                                    :: ehfx1, ehfx2
     REAL(KIND=dp)                                    :: xfraction
     TYPE(semi_empirical_si_type), POINTER            :: gem_integral_storage
     TYPE(semi_empirical_si_type), POINTER            :: gemdiag_integral_storage
     TYPE(semi_empirical_si_type), POINTER            :: ggg_integral_storage
     REAL(KIND=dp)                                    :: screening
     REAL(KIND=dp)                                    :: eps_opt
     INTEGER                                          :: max_iter
  END TYPE ri_environment_type

  TYPE ri_environment_p_type
     TYPE(ri_environment_type),POINTER :: ri_env
  END TYPE ri_environment_p_type

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_environment_types'
  INTEGER, SAVE, PRIVATE :: last_ri_env_id=0

  PUBLIC :: ri_environment_type,&
            ri_vector_type,&
            set_ri_env, ri_env_create,&
            get_ri_env, ri_env_release,&
            ri_vector_create, ri_vector_release,&
            ri_vector_add, ri_vector_copy, ri_vector_set,&
            ri_vector_dot, ri_vector_norm, ri_vector_scale,&
            ri_vector_sync

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an ri_env
!> \param ri_env the ri_environment you want to create
!>
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE ri_env_create ( ri_env, error )

    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( ri_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_ri_env_id=last_ri_env_id+1
    ri_env%id_nr=last_ri_env_id
    ri_env%ref_count=1
    ri_env%in_use=0

    ri_env%ehfx1=0._dp
    ri_env%ehfx2=0._dp
    ri_env%xfraction=1._dp

    ri_env%screening=0._dp
    ri_env%eps_opt=0._dp
    ri_env%max_iter=0

    NULLIFY(ri_env%gem_integral_storage)
    NULLIFY(ri_env%gemdiag_integral_storage)
    NULLIFY(ri_env%ggg_integral_storage)

  END SUBROUTINE ri_env_create

! *****************************************************************************
!> \brief provides a method for attaching various structures to an ri_env
!> \param ri_env the ri_environment you want to change
!> \param coeff ...
!> \param rhs ...
!> \param gem_integral_storage ...
!> \param ggg_integral_storage ...
!> \param gemdiag_integral_storage ...
!> \param xfraction ...
!> \param screening ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE set_ri_env ( ri_env, coeff, rhs, gem_integral_storage, ggg_integral_storage, &
                          gemdiag_integral_storage, xfraction, screening, error )

    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), OPTIONAL, POINTER  :: coeff, rhs
    TYPE(semi_empirical_si_type), OPTIONAL, &
      POINTER                                :: gem_integral_storage, &
                                                ggg_integral_storage, &
                                                gemdiag_integral_storage
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: xfraction, screening
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    IF ( PRESENT(coeff) ) ri_env%coeff => coeff
    IF ( PRESENT(rhs) ) ri_env%rhs => rhs
    IF ( PRESENT(gem_integral_storage) ) THEN
       CALL semi_empirical_si_release(ri_env%gem_integral_storage,error=error)
       ri_env%gem_integral_storage => gem_integral_storage
    END IF
    IF ( PRESENT(ggg_integral_storage) ) THEN
       CALL semi_empirical_si_release(ri_env%ggg_integral_storage,error=error)
       ri_env%ggg_integral_storage => ggg_integral_storage
    END IF
    IF ( PRESENT(gemdiag_integral_storage) ) THEN
       CALL semi_empirical_si_release(ri_env%gemdiag_integral_storage,error=error)
       ri_env%gemdiag_integral_storage => gemdiag_integral_storage
    END IF
    IF ( PRESENT(screening) ) ri_env%screening = screening
    IF ( PRESENT(xfraction) ) ri_env%xfraction = xfraction

  END SUBROUTINE set_ri_env

! *****************************************************************************
!> \brief provides a method for getting the various structures attached
!>      to an ri_env
!> \param ri_env the ri_environment you want to get information on
!> \param coeff ...
!> \param rhs ...
!> \param gem_integral_storage ...
!> \param ggg_integral_storage ...
!> \param gemdiag_integral_storage ...
!> \param screening ...
!> \param xfraction ...
!> \param eps_opt ...
!> \param max_iter ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE get_ri_env ( ri_env, coeff, rhs, gem_integral_storage, ggg_integral_storage, &
                          gemdiag_integral_storage, screening, xfraction, eps_opt, max_iter, error )

    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), OPTIONAL, POINTER  :: coeff, rhs
    TYPE(semi_empirical_si_type), OPTIONAL, &
      POINTER                                :: gem_integral_storage, &
                                                ggg_integral_storage, &
                                                gemdiag_integral_storage
    REAL(KIND=dp), INTENT(inout), OPTIONAL   :: screening, xfraction, eps_opt
    INTEGER, INTENT(inout), OPTIONAL         :: max_iter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    IF ( PRESENT(coeff) ) coeff => ri_env%coeff
    IF ( PRESENT(rhs) ) rhs => ri_env%rhs
    IF ( PRESENT(gem_integral_storage) ) gem_integral_storage => ri_env%gem_integral_storage
    IF ( PRESENT(ggg_integral_storage) ) ggg_integral_storage => ri_env%ggg_integral_storage
    IF ( PRESENT(gemdiag_integral_storage) ) gemdiag_integral_storage => ri_env%gemdiag_integral_storage
    IF ( PRESENT(screening) ) screening = ri_env%screening
    IF ( PRESENT(xfraction) ) xfraction = ri_env%xfraction
    IF ( PRESENT(eps_opt) ) eps_opt = ri_env%eps_opt
    IF ( PRESENT(max_iter) ) max_iter = ri_env%max_iter

  END SUBROUTINE get_ri_env

! *****************************************************************************
!> \brief retains the given ri env
!> \param ri_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!> \author JGH
! *****************************************************************************
SUBROUTINE ri_env_retain(ri_env, error)
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(ri_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(ri_env%ref_count>0,cp_failure_level,routineP,error)
     ri_env%ref_count=ri_env%ref_count+1
  END IF

END SUBROUTINE ri_env_retain

! *****************************************************************************
!> \brief releases the given mc env
!> \param ri_env the mc environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!> \author JGH
! *****************************************************************************
SUBROUTINE ri_env_release(ri_env, error)
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(ri_env)) THEN
     CPPreconditionNoFail(ri_env%ref_count>0,cp_failure_level,routineP,error)
     ri_env%ref_count=ri_env%ref_count-1
     IF (ri_env%ref_count==0) THEN
        ri_env%ref_count=1
        ri_env%ref_count=0
        CALL ri_vector_release(ri_env%coeff, error)
        CALL ri_vector_release(ri_env%rhs, error)
        IF (ASSOCIATED(ri_env%gem_integral_storage)) THEN
           CALL semi_empirical_si_release(ri_env%gem_integral_storage, error)
        END IF
        IF (ASSOCIATED(ri_env%ggg_integral_storage)) THEN
           CALL semi_empirical_si_release(ri_env%ggg_integral_storage, error)
        END IF
        IF (ASSOCIATED(ri_env%gemdiag_integral_storage)) THEN
           CALL semi_empirical_si_release(ri_env%gemdiag_integral_storage, error)
        END IF
        DEALLOCATE(ri_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(ri_env)
END SUBROUTINE ri_env_release

! *****************************************************************************
!> \brief ...
!> \param ri_vec ...
!> \param nkind ...
!> \param natoms ...
!> \param nbasis ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_create(ri_vec, nkind, natoms, nbasis, error)
    TYPE(ri_vector_type), POINTER            :: ri_vec
    INTEGER, INTENT(IN)                      :: nkind
    INTEGER, DIMENSION(:), INTENT(IN)        :: natoms, nbasis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, tl

   CALL ri_vector_release(ri_vec, error)

   ALLOCATE(ri_vec,stat=istat)
   CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
   ALLOCATE(ri_vec%vector(nkind),stat=istat)
   CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

   tl=0
   DO ikind=1,nkind
     ri_vec%vector(ikind)%natom = natoms(ikind)
     ri_vec%vector(ikind)%nbasis = nbasis(ikind)
     tl=tl+natoms(ikind)*nbasis(ikind)
     ALLOCATE(ri_vec%vector(ikind)%vmat(nbasis(ikind),natoms(ikind)),stat=istat)
     CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
     ri_vec%vector(ikind)%vmat=0._dp
   END DO
   ri_vec%total_length = tl

END SUBROUTINE ri_vector_create

! *****************************************************************************
!> \brief ...
!> \param ri_vec ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_release(ri_vec, error)
    TYPE(ri_vector_type), POINTER            :: ri_vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat

   IF(ASSOCIATED(ri_vec)) THEN
     IF(ASSOCIATED(ri_vec%vector)) THEN
       DO ikind=1,SIZE(ri_vec%vector)
          DEALLOCATE(ri_vec%vector(ikind)%vmat,stat=istat)
          CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
       END DO
       DEALLOCATE(ri_vec%vector,stat=istat)
       CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
       NULLIFY(ri_vec%vector)
     END IF
     DEALLOCATE(ri_vec,stat=istat)
     CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
   END IF

END SUBROUTINE ri_vector_release
! *****************************************************************************
!> \brief ...
!> \param ri_vec ...
!> \param value ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_set(ri_vec, value, error)
    TYPE(ri_vector_type)                     :: ri_vec
    REAL(KIND=dp), INTENT(IN)                :: value
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   DO i=1,SIZE(ri_vec%vector)
     ri_vec%vector(i)%vmat(:,:) = value
   END DO

END SUBROUTINE ri_vector_set

! *****************************************************************************
!> \brief ...
!> \param alpha ...
!> \param vecx ...
!> \param vecy ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_add(alpha, vecx, vecy, error)
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(ri_vector_type)                     :: vecx, vecy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   DO i=1,SIZE(vecx%vector)
     vecy%vector(i)%vmat(:,:) = vecy%vector(i)%vmat(:,:) + alpha*vecx%vector(i)%vmat(:,:)
   END DO

END SUBROUTINE ri_vector_add

! *****************************************************************************
!> \brief ...
!> \param vec_in ...
!> \param vec_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_copy(vec_in, vec_out, error)
    TYPE(ri_vector_type)                     :: vec_in, vec_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   DO i=1,SIZE(vec_in%vector)
     vec_out%vector(i)%vmat(:,:) = vec_in%vector(i)%vmat(:,:)
   END DO

END SUBROUTINE ri_vector_copy

! *****************************************************************************
!> \brief ...
!> \param scale ...
!> \param vec ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_scale(scale, vec, error)
    REAL(KIND=dp), INTENT(IN)                :: scale
    TYPE(ri_vector_type)                     :: vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   DO i=1,SIZE(vec%vector)
     vec%vector(i)%vmat(:,:) = scale*vec%vector(i)%vmat(:,:)
   END DO

END SUBROUTINE ri_vector_scale

! *****************************************************************************
!> \brief ...
!> \param RESULT ...
!> \param vec1 ...
!> \param vec2 ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_dot(RESULT, vec1, vec2, error)
    REAL(KIND=dp), INTENT(OUT)               :: RESULT
    TYPE(ri_vector_type)                     :: vec1, vec2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_dot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   RESULT = 0._dp
   DO i=1,SIZE(vec1%vector)
     RESULT = RESULT + SUM(vec1%vector(i)%vmat(:,:)*vec2%vector(i)%vmat(:,:))
   END DO

END SUBROUTINE ri_vector_dot

! *****************************************************************************
!> \brief ...
!> \param RESULT ...
!> \param vec ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_norm(RESULT, vec, error)
    REAL(KIND=dp), INTENT(OUT)               :: RESULT
    TYPE(ri_vector_type)                     :: vec
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   RESULT = 0._dp
   DO i=1,SIZE(vec%vector)
     RESULT = RESULT + SUM(vec%vector(i)%vmat(:,:)*vec%vector(i)%vmat(:,:))
   END DO
   RESULT = SQRT ( RESULT/REAL(vec%total_length,KIND=dp) )

END SUBROUTINE ri_vector_norm


! *****************************************************************************
!> \brief ...
!> \param vec ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE ri_vector_sync(vec, para_env, error)
    TYPE(ri_vector_type)                     :: vec
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ri_vector_sync', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

   DO i=1,SIZE(vec%vector)
     CALL mp_sum(vec%vector(i)%vmat(:,:),para_env%group)
   END DO

END SUBROUTINE ri_vector_sync
! *****************************************************************************

END MODULE ri_environment_types

