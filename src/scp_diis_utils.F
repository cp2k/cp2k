!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE scp_diis_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE linear_systems,                  ONLY: solve_system
  USE message_passing,                 ONLY: mp_max
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE termination,                     ONLY: stop_memory
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_diis_utils'
! *****************************************************************************
  TYPE diis_info_type
    INTEGER :: nhist, ihist
    REAL (KIND=dp), DIMENSION(:,:), POINTER :: e_hist
    REAL (KIND=dp), DIMENSION(:,:), POINTER :: c_hist
    REAL (KIND=dp), DIMENSION(:), POINTER :: hess
  END TYPE diis_info_type


CONTAINS

! *****************************************************************************
!> \brief ...
!> \param diis_info ...
!> \param nhist ...
!> \param aux_coeff_set ...
! *****************************************************************************
  SUBROUTINE diis_info_create(diis_info,nhist,aux_coeff_set)

    TYPE(diis_info_type), POINTER            :: diis_info
    INTEGER, INTENT(in)                      :: nhist
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set

    CHARACTER(len=*), PARAMETER :: routineN = 'diis_info_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, ncoef, ncoef_of_kind, &
                                                stat
    LOGICAL                                  :: failure
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cp_error_type)                      :: error

    CPPrecondition(ASSOCIATED(aux_coeff_set),cp_failure_level,routineP,error,failure)
    ALLOCATE (diis_info,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "diis_info",0)
    diis_info%nhist=nhist
    diis_info%ihist=0
    ncoef=0
    NULLIFY(local_coeffs,diis_info%e_hist,diis_info%c_hist, diis_info%hess)
    DO ikind=1,SIZE(aux_coeff_set%coeffs_of_kind)
      local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF(ASSOCIATED(local_coeffs))THEN
        CALL get_aux_coeff(coeffs=local_coeffs,&
                           ncoef_of_kind=ncoef_of_kind,error=error)
        ncoef=ncoef + ncoef_of_kind
      END IF
    END DO
    IF (.not.ASSOCIATED(diis_info%c_hist)) THEN
      ALLOCATE (diis_info% c_hist(nhist,ncoef),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                      "diis_info%c_hist",dp_size*ncoef*nhist)
    END IF
    IF (.not.ASSOCIATED(diis_info%e_hist)) THEN
      ALLOCATE (diis_info% e_hist(nhist,ncoef),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                      "diis_info%e_hist",dp_size*ncoef*nhist)
    END IF
    IF (.not.ASSOCIATED(diis_info%hess)) THEN
      ALLOCATE (diis_info%hess(ncoef),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                      "diis_info%hess",dp_size*ncoef)
    END IF
    diis_info%c_hist(:,:)=0.0_dp
    diis_info%e_hist(:,:)=0.0_dp
    diis_info%hess(:)=1.0_dp

  END SUBROUTINE diis_info_create

! *****************************************************************************
!> \brief ...
!> \param diis_info ...
! *****************************************************************************
  SUBROUTINE diis_info_destroy(diis_info)

    TYPE(diis_info_type), POINTER            :: diis_info

    CHARACTER(len=*), PARAMETER :: routineN = 'diis_info_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(diis_info)) THEN
      IF (ASSOCIATED(diis_info%c_hist)) THEN
        DEALLOCATE (diis_info% c_hist,STAT=stat)
        IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"diis_info%c_hist")
      END IF
      IF (ASSOCIATED(diis_info%e_hist)) THEN
        DEALLOCATE (diis_info%e_hist,STAT=stat)
        IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"diis_info%e_hist")
      END IF
      IF (ASSOCIATED(diis_info%hess)) THEN
        DEALLOCATE (diis_info%hess,STAT=stat)
        IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"diis_info%hess")
      END IF
      DEALLOCATE (diis_info,STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"diis_info")
    END IF
    NULLIFY(diis_info)

  END SUBROUTINE diis_info_destroy

! *****************************************************************************
!> \brief ...
!> \param scp_env ...
!> \param atomic_kind_set ...
!> \param diis_info ...
!> \param maxerr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE do_diis_step ( scp_env, atomic_kind_set, diis_info, maxerr, error)

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(diis_info_type), POINTER            :: diis_info
    REAL(KIND=dp), INTENT(INOUT)             :: maxerr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_diis_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, icoef, icoef_atom, ikind, &
                                                j, n_els, ncoef, ncoef_atom, &
                                                ncoef_of_kind, nhist, nkind, &
                                                stat
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: c_forces, cg, dcg, eg
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c, fc
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs

  NULLIFY(aux_coeff_set,c,fc,local_coeffs)
  CALL get_scp_env(scp_env=scp_env, &
                  aux_coeff_set=aux_coeff_set,error=error)

  nkind = SIZE (atomic_kind_set)
  nhist = diis_info%nhist
  ncoef = 0
  icoef = 0
  maxerr = 0.0_dp

  DO ikind=1,nkind
     local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
     IF(ASSOCIATED(local_coeffs))THEN
      CALL get_aux_coeff(coeffs=local_coeffs, ncoef_of_kind=ncoef_of_kind,error=error)
       ncoef=ncoef + ncoef_of_kind
     END IF
  END DO
 IF (ncoef/=0) THEN
! allocating the arrays
  ALLOCATE (cg(ncoef),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                  "cg",dp_size*ncoef)
  ALLOCATE (eg(ncoef),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                  "eg",dp_size*ncoef)
  ALLOCATE (c_forces(ncoef),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                  "c_forces",dp_size*ncoef)
  ALLOCATE (dcg(nhist+1),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                  "dcg",dp_size*(nhist+1))

! increment the counters for the histories and the iterations
  diis_info%ihist=diis_info%ihist+1

! forming the trial vector from the natural force
  DO ikind= 1, nkind
     local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
     IF(ASSOCIATED(local_coeffs))THEN
       CALL get_aux_coeff(coeffs=local_coeffs,c=c,fc=fc,&
                          n_els=n_els,ncoef_atom=ncoef_atom,error=error)
       DO i=1,n_els
          DO icoef_atom=1,ncoef_atom
             icoef = icoef + 1
             cg (icoef) = c(i,icoef_atom)
             c_forces(icoef)=fc(i,icoef_atom)
             WRITE ( * , * ) 'FCOEF', ikind, i, fc ( i, icoef_atom )
             c(i,icoef_atom)=0.0_dp
          END DO
       END DO
     END IF
  END DO

! store the histories
  CALL store_hist(c_forces, cg, diis_info%e_hist, diis_info%c_hist,&
                   diis_info%hess, nhist, diis_info%ihist)

! solve for the coefficients
  CALL solve(dcg,diis_info%e_hist,diis_info%ihist,ncoef)

! the new coefficients!
  eg(:)=0.0_dp
  DO j=1, diis_info%ihist
     icoef = 0
     DO ikind= 1, nkind
        local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
        IF(ASSOCIATED(local_coeffs))THEN
          CALL get_aux_coeff(coeffs=local_coeffs,c=c,&
                             n_els=n_els,ncoef_atom=ncoef_atom,error=error)
          DO i=1,n_els
             DO icoef_atom = 1, ncoef_atom
                icoef=icoef+1
                c(i,icoef_atom)=c(i,icoef_atom) + dcg(j)*diis_info%c_hist(j,icoef)
                eg (icoef) = eg(icoef) + dcg(j) * diis_info%e_hist(j,icoef)
             END DO
          END DO
        END IF
     ENDDO
  ENDDO
  maxerr=MAXVAL(ABS(eg))
 ELSE
  maxerr=0
 END IF
 CALL mp_max(maxerr,scp_env%para_env%group)
! estimate new coefficents
 IF (ncoef/=0) THEN
  icoef=0
  DO ikind= 1, nkind
    local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(ASSOCIATED(local_coeffs))THEN
      CALL get_aux_coeff(coeffs=local_coeffs,c=c,&
                          n_els=n_els,ncoef_atom=ncoef_atom,error=error)
      DO i=1,n_els
         DO icoef_atom = 1, ncoef_atom
            icoef=icoef+1
            c(i,icoef_atom)=c(i,icoef_atom) - eg(icoef)*diis_info%hess(icoef)
            WRITE ( * , * ) 'COEF', ikind, i, c ( i, icoef_atom )
         END DO
      END DO
    END IF
  ENDDO
! deallocating the arrays
  DEALLOCATE (cg,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cg")
  DEALLOCATE (c_forces,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"c_forces")
  DEALLOCATE (eg,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eg")
  DEALLOCATE (dcg,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dcg")
  END IF

  END SUBROUTINE do_diis_step

! *****************************************************************************
!> \brief ...
!> \param c_forces ...
!> \param cg ...
!> \param e_hist ...
!> \param c_hist ...
!> \param hess ...
!> \param nhist ...
!> \param ihist ...
! *****************************************************************************
  SUBROUTINE STORE_HIST(c_forces,cg,e_hist,c_hist,hess,nhist,ihist)

    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: c_forces
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: cg
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: e_hist, c_hist
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: hess
    INTEGER, INTENT(IN)                      :: nhist
    INTEGER, INTENT(INOUT)                   :: ihist

    INTEGER                                  :: icoef, j, ncoef

  ncoef = SIZE(c_forces)

! constructing the error vector and the coefficient vector over histories
  IF (nhist==1) THEN
    ihist=nhist
    GOTO 2000
  ELSE
    IF (ihist>nhist) THEN
      ihist=nhist
      DO icoef=1,ncoef
        DO j=2, nhist
          e_hist(j-1,icoef)=e_hist(j,icoef)
          c_hist(j-1,icoef)=c_hist(j,icoef)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
2000  CONTINUE
  DO icoef = 1, ncoef
    c_hist (ihist,icoef) = cg (icoef)
    e_hist (ihist,icoef) =  c_forces (icoef) / hess(icoef)
  ENDDO

  END SUBROUTINE store_hist

! *****************************************************************************
!> \brief ...
!> \param dcg ...
!> \param e_hist ...
!> \param ihist ...
!> \param ncoef ...
! *****************************************************************************
  SUBROUTINE SOLVE(dcg,e_hist,ihist,ncoef)
    REAL(KIND=dp), INTENT(INOUT)             :: dcg(:)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e_hist
    INTEGER, INTENT(IN)                      :: IHIST, ncoef

    CHARACTER(len=*), PARAMETER :: routineN = 'SOLVE', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icoef, j, k, mysize, stat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: bcg, eigenvec

! constructing the solution vector

  DO j=1, ihist
   dcg(j)=0.0_dp
  ENDDO
  dcg(ihist+1)=-1.0_dp

! constructing the overlap matrix
  ALLOCATE (bcg(ihist+1,ihist+1),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                  "bcg",dp_size*(ihist+1)**2)
  ALLOCATE (eigenvec(ihist+1,1),STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                  "eigenvec",dp_size*(ihist+1))

  bcg=0.0_dp
  DO j=1,ihist
    DO k=1,ihist
      DO icoef = 1, ncoef
        bcg(j,k) = bcg(j,k) + e_hist(j,icoef) * e_hist(k,icoef)
      ENDDO
    ENDDO
  ENDDO
  DO j=1, ihist
    bcg(ihist+1,j) = -1.0_dp
    bcg(j,ihist+1) = -1.0_dp
  ENDDO
  DO j=1, ihist + 1
    eigenvec(j,1) = dcg(j)
  END DO

  bcg (ihist+1,ihist+1) = 0.0_dp
  mysize = ihist + 1

! get the solution vector

  CALL solve_system ( bcg, mysize, eigenvec)

  DO j=1, ihist
     dcg(j) =  eigenvec(j,1)
  END DO

  DEALLOCATE (bcg,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"bcg")
  DEALLOCATE (eigenvec,STAT=stat)
  IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eigenvec")

  END SUBROUTINE solve

END MODULE scp_diis_utils

