!------------------------------------------------------------------------------!
  MODULE pol_diis
!------------------------------------------------------------------------------!
!
  USE kinds,                           ONLY: dbl
  USE linear_systems,                  ONLY: solve_system
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
!
  PUBLIC :: diis
  PRIVATE
!      
!------------------------------------------------------------------------------!
!
  CONTAINS
!
!------------------------------------------------------------------------------!
  SUBROUTINE DIIS ( c_pos, c_force, dtw, fconv, hess, &
                   e_hist, c_hist, nhist, ihist)
!-----
!  c_pos: basis functions coefficents
!  c_force: forces on the basis functions coefficents
!  dtw: is the time-step, dt
!  fconv: maximum gradient       
!  hess: is the Hessian  
!  e_hist: is a history vector
!  c_hist: is the history vector of the coefficients 
!  nhist: number of histories to be stored
!  ihist: history counter
!--------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: c_pos, c_force
    REAL(dbl), INTENT(IN)                    :: dtw
    REAL(dbl), INTENT(INOUT)                 :: fconv
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: hess
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: e_hist, c_hist
    INTEGER, INTENT(IN)                      :: nhist
    INTEGER, INTENT(INOUT)                   :: ihist

    INTEGER                                  :: icoef, j, ncoef
    REAL(dbl), ALLOCATABLE                   :: dcg(:)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: cg, eg

! allocating the arrays

   ncoef = size(c_pos)
   ALLOCATE(cg(ncoef)) 
   ALLOCATE(eg(ncoef))
   ALLOCATE(dcg(nhist+1))

! increment the counters for the histories and the iterations
  ihist=ihist+1

! forming the trial vector from the natural force 
  DO icoef = 1, ncoef
    cg (icoef) = c_pos (icoef) + dtw/hess(icoef) * c_force (icoef)
  END DO

! store the histories
  CALL store_hist(c_force, dtw, cg, e_hist, c_hist, hess, nhist, ihist)

! solve for the coefficients
  CALL solve(dcg,e_hist,ihist,ncoef)

! the new coefficients!
  DO icoef = 1, ncoef
    c_pos (icoef) = 0._dbl
  END DO
  eg(:)=0._dbl
  DO j=1, ihist
    DO icoef = 1, ncoef
      c_pos (icoef) = c_pos (icoef) + dcg(j) * c_hist(j,icoef)
      eg (icoef) = eg(icoef) + dcg(j) * e_hist(j,icoef)
    ENDDO
  ENDDO
  fconv=maxval(abs(eg))
! deallocating the arrays
  DEALLOCATE(cg)
  DEALLOCATE(eg)
  DEALLOCATE(dcg)
  return
  END SUBROUTINE diis
!------------------------------------------------------------------------------!
  SUBROUTINE STORE_HIST(c_forces,dtw,cg,e_hist,c_hist,hess,nhist,ihist)
!------------------------------------------------------------------------------!

    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: c_forces
    REAL(dbl), INTENT(IN)                    :: dtw
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: cg
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: e_hist, c_hist
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: hess
    INTEGER, INTENT(IN)                      :: nhist
    INTEGER, INTENT(INOUT)                   :: ihist

    INTEGER                                  :: icoef, j, ncoef

  ncoef = size(c_forces)

! constructing the error vector and the coefficient vector over histories
  IF (nhist==1) THEN
    ihist=nhist
    goto 2000
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
2000  continue
  DO icoef = 1, ncoef
    c_hist (ihist,icoef) = cg (icoef)
    e_hist (ihist,icoef) = dtw * c_forces (icoef) / hess(icoef)
  ENDDO
  return
  END SUBROUTINE store_hist
!------------------------------------------------------------------------------!
  SUBROUTINE SOLVE(dcg,e_hist,ihist,ncoef)
!------------------------------------------------------------------------------!
    REAL(dbl), INTENT(OUT)                   :: dcg(:)
    REAL(dbl), DIMENSION(:, :), intent(IN)   :: e_hist
    INTEGER, INTENT(IN)                      :: IHIST, ncoef

    INTEGER                                  :: icoef, ios, j, k, mysize
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: bcg, eigenvec

! constructing the solution vector

  do j=1, ihist
   dcg(j)=0._dbl
  enddo
  dcg(ihist+1)=-1._dbl

! constructing the overlap matrix
  ALLOCATE (bcg(ihist+1,ihist+1), STAT = ios)
  IF ( ios /= 0 )  CALL stop_memory ( 'pol_diis', 'bcg', ihist +1 )
  ALLOCATE (eigenvec(ihist+1,1), STAT = ios)
  IF ( ios /= 0 )  CALL stop_memory ( 'pol_diis', 'eigenvec', ihist +1 )

  bcg=0._dbl
  do j=1,ihist
    do k=1,ihist
      do icoef = 1, ncoef
        bcg(j,k) = bcg(j,k) + e_hist(j,icoef) * e_hist(k,icoef)
      enddo
    enddo
  enddo
  do j=1, ihist
    bcg(ihist+1,j) = -1._dbl
    bcg(j,ihist+1) = -1._dbl
  enddo
  do j=1, ihist + 1
    eigenvec(j,1) = dcg(j)
  end do

  bcg (ihist+1,ihist+1) = 0._dbl
  mysize = ihist + 1

! get the solution vector

  call solve_system ( bcg, mysize, eigenvec)
 
  do j=1, ihist
     dcg(j) =  eigenvec(j,1) 
  end do

  DEALLOCATE(bcg)
  IF ( ios /= 0 )  CALL stop_memory ( 'pol_diis', 'deall bcg' )
  DEALLOCATE(eigenvec)
  IF ( ios /= 0 )  CALL stop_memory ( 'pol_diis', 'deall eigenvec' )
  RETURN
  END SUBROUTINE solve
!------------------------------------------------------------------------------!
  END MODULE pol_diis
!------------------------------------------------------------------------------!
