!-----------------------------------------------------------------------------!
!   FES: a fast and general program to map metadynamics on grids              !
!   Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009 Teodoro Laino       !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Program to Map on grid the hills spawned during a metadynamics run
!> \author Teodoro Laino [tlaino] - 06.2009
!> \par History
!>     03.2006 created [tlaino]
!>     teodoro.laino .at. gmail.com 
!>     11.2007 - tlaino (University of Zurich): Periodic COLVAR - cleaning.
!>
!> \par Note
!>     Please report any bug to the author
! *****************************************************************************
MODULE graph_methods

  USE graph_utils,                     ONLY: derivative,&
                                             point_no_pbc,&
                                             point_pbc
  USE kinds,                           ONLY: dp

  IMPLICIT NONE
  PUBLIC :: fes_compute_low,&
            fes_write,&
            fes_only_write,&
            fes_min,&
            fes_path

CONTAINS
! *****************************************************************************
!> \brief Efficiently map the gaussians on the grid
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE fes_compute_low(idim, nn, fes, gauss, ind, ind0, nfes, ndim,&
       ngauss, ngrid, iperd)
    INTEGER, INTENT(in)                      :: idim
    INTEGER, DIMENSION(:)                    :: nn
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gauss
    INTEGER, DIMENSION(:)                    :: ind, ind0
    INTEGER, INTENT(in)                      :: nfes, ndim, ngauss
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    INTEGER, DIMENSION(:)                    :: iperd

    INTEGER                                  :: i, j, k, pnt
    INTEGER, DIMENSION(:), POINTER           :: ll, pos
    REAL(KIND=dp)                            :: prod

    ALLOCATE(pos(ndim),ll(ndim))
    pos=ind
    k=nn(idim)

    DO i = -k, k
       pos(idim) = ind(idim)+i
       IF (iperd(idim) == 0 ) THEN
          IF (pos(idim) .GT. ngrid(idim)) CYCLE
          IF (pos(idim) .LT. 1) CYCLE
       END IF
       IF(idim/=1) THEN 
          CALL fes_compute_low(idim-1, nn, fes, gauss, pos, ind0, nfes, ndim, ngauss, ngrid, iperd)
       ELSE
          pnt  = point_pbc(pos,iperd,ngrid,ndim)
          prod = 1.0_dp
          DO j=1,ndim
             ll(j) = pos(j) - ind0(j)
             prod  = prod*gauss(ll(j),j)
          END DO
          fes(pnt)=fes(pnt)+prod
       END IF
    END DO
    DEALLOCATE(pos,ll)

  END SUBROUTINE fes_compute_low

! *****************************************************************************
!> \brief Writes the FES on the file
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw)
    INTEGER, INTENT(IN)                      :: idim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, DIMENSION(:), POINTER           :: pos
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid, x0
    INTEGER, INTENT(IN)                      :: ndw

    INTEGER                                  :: dimval, i, id, pnt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xx

    ALLOCATE(xx(ndim))
    xx = x0
    DO i = 1,ngrid(idim)
       pos(idim)=i
       IF(idim/=ndim-ndw+1) THEN
          CALL fes_write(idim-1, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw)
       ELSE        
          pnt=point_no_pbc(pos,ngrid,ndim)
          xx = x0 + dp_grid * ( pos - 1 )
          dimval = PRODUCT(ngrid(1:ndim-ndw)) 
          WRITE(123,'(10f20.10)')(xx(id),id=ndim,ndim-ndw+1,-1),MINVAL(-fes(pnt:pnt+dimval-1))
       END IF
    END DO
    DEALLOCATE(xx)

  END SUBROUTINE fes_write

! *****************************************************************************
!> \brief Writes the FES on the file when stride is requested
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE fes_only_write(idim, fes,  pos,  ndim, ngrid, ndw)
    INTEGER, INTENT(IN)                      :: idim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, DIMENSION(:), POINTER           :: pos
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    INTEGER, INTENT(IN)                      :: ndw

    INTEGER                                  :: dimval, i, pnt

    DO i = 1,ngrid(idim)
       pos(idim)=i
       IF(idim/=ndim-ndw+1) THEN
          CALL fes_only_write(idim-1, fes,  pos,  ndim, ngrid, ndw)
       ELSE        
          pnt=point_no_pbc(pos,ngrid,ndim)
          dimval = PRODUCT(ngrid(1:ndim-ndw))
          WRITE(123,'(1f12.5)')MINVAL(-fes(pnt:pnt+dimval-1))
       END IF
    END DO

  END SUBROUTINE fes_only_write

! *****************************************************************************
!> \brief Finds minima of the FES
!> \par History
!>      06.2009 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fes_min(fes, ndim, ngrid, dp_grid, x0, ndw)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid, x0
    INTEGER, INTENT(IN)                      :: ndw

    INTEGER                                  :: i, id, j, k, nacc, ntrials, &
                                                pnt
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: history
    INTEGER, DIMENSION(:), POINTER           :: pos, pos0
    INTEGER, DIMENSION(ndim)                 :: Dpos, ntrust
    LOGICAL                                  :: do_save
    REAL(KIND=dp)                            :: norm_dx, resto
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dx, rnd, xx

    ntrust  = ngrid/10
    ntrials = PRODUCT (ngrid)
    WRITE(6,'(A,10I6)',ADVANCE="no")"FES| Trust hyper-radius ",ntrust
    WRITE(6,'(A,10F12.6)')" which is equivalent to: ",ntrust*dp_grid

    ALLOCATE(xx(ndim),dx(ndim),pos0(ndim),rnd(ndim),pos(ndim))
    ALLOCATE(history(ndim,ntrials))
    history = 0
    nacc    = 0
    DO j = 1, ntrials
       ! Loop over all points
       pnt = j
       DO k = ndim, 2, -1
          pos0(k) = pnt/PRODUCT(ngrid(1:k-1))
          resto = MOD(pnt,PRODUCT(ngrid(1:k-1)))
          IF (resto/=0) THEN 
             pnt = pnt -pos0(k)*PRODUCT(ngrid(1:k-1))
             pos0(k)=pos0(k)+1
          ELSE
             pnt = PRODUCT(ngrid(1:k-1))
          END IF
       END DO
       pos0(1) = pnt

       IF (ANY(pos0<ntrust))       CYCLE
       IF (ANY(pos0>ngrid-ntrust)) CYCLE

       ! Evaluate position and derivative
       pos = pos0
       xx  = x0 + dp_grid * ( pos - 1 )
       dx = derivative (fes, pos, ndim, ngrid, dp_grid)

       ! Integrate till derivative is small enough..
       DO i = 1, 100
          pnt = point_no_pbc(pos,ngrid,ndim)
          !WRITE(10+j,'(10f20.10)')(xx(id),id=ndim,1,-1),-fes(pnt)  

          norm_dx = SQRT(DOT_PRODUCT(dx,dx))
          IF (norm_dx == 0.0_dp) EXIT ! It is in a really flat region
          xx  = xx - MIN(0.1_dp,norm_dx)*dx/norm_dx
          ! Re-evaluating pos
          pos = CEILING((xx-x0)/dp_grid)+1

          ! Incremental pos
          dx  = derivative (fes, pos, ndim, ngrid, dp_grid)
       END DO

       ! Compare with the available minima and if they are the same skip
       ! saving this position..
       do_save = .TRUE.
       DO i = 1, nacc
          Dpos = pos - history(:,i)
          norm_dx  = DOT_PRODUCT(Dpos,Dpos)
          IF ((SQRT(REAL(norm_dx,KIND=dp))<=MAXVAL(ntrust)).OR.(fes(pnt)<1.0E-3_dp)) THEN 
             do_save = .FALSE.
             EXIT
          END IF
       END DO
       IF (do_save) THEN
          pnt = point_no_pbc(pos,ngrid,ndim)
          xx  = x0+dp_grid*(pos-1)
          WRITE(6,'(A,5F12.6)',ADVANCE="NO") "FES| Minimum found (",(xx(id),id=ndim,ndim-ndw+1,-1)
          WRITE(6,'(A,F12.6,A)')" ). FES value = ",-fes(pnt)," Hartree."
          nacc = nacc + 1
          history (:,nacc) = pos
       END IF
    END DO
    WRITE(6,'(A,I6,A)') "FES| Number of Minimum found: ",nacc,"."

    DEALLOCATE(xx,dx,pos0,rnd,pos)
    DEALLOCATE(history)

  END SUBROUTINE fes_min

! *****************************************************************************
!> \brief Finds path between two points (a) and (b)
!> \par History
!>      06.2009 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fes_path(fes, ndim, ngrid, dp_grid, x0, ndw, minima)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid, x0
    INTEGER, INTENT(IN)                      :: ndw
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: minima

  END SUBROUTINE fes_path

END MODULE graph_methods
