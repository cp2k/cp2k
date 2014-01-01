!-----------------------------------------------------------------------------!
!   FES: a fast and general program to map metadynamics on grids              !
!   Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013 !
!                 Teodoro Laino                                               !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Program to Map on grid the hills spawned during a metadynamics run
!> \author Teodoro Laino [tlaino] - 06.2009
!> \par History
!>     03.2006 created [tlaino]
!>     teodoro.laino .at. gmail.com 
!>     11.2007 - tlaino (University of Zurich): Periodic COLVAR - cleaning.
!>     12.2010 - teodoro.laino@gmail.com: addition of the MEP for FES
!>
!> \par Note
!>     Please report any bug to the author
! *****************************************************************************
MODULE graph_methods

  USE graph_utils,                     ONLY: derivative,&
                                             mep_input_data_type,&
                                             pbc,&
                                             point_no_pbc,&
                                             point_pbc
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE physcon,                         ONLY: bohr, init_physcon
  USE string_utilities,                ONLY: uppercase
  USE periodic_table,                  ONLY: ptable, nelem
  USE graph_utils,                     ONLY: get_val_res

  IMPLICIT NONE
  PUBLIC :: fes_compute_low,&
            fes_write,&
            fes_only_write,&
            fes_min,&
            fes_path,&
            fes_cube_write

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
  RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid,&
       & dp_grid, x0, ndw, l_fes_int, array)
    INTEGER, INTENT(IN)                      :: idim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, DIMENSION(:), POINTER           :: pos
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid, x0
    INTEGER, INTENT(IN)                      :: ndw
    LOGICAL, INTENT(IN)                      :: l_fes_int
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: array

    INTEGER                                  :: dimval, i, id, is,&
         & pnt, ind, it, itt, np
    REAL(KIND=dp)                            :: dvol, sum_fes
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xx

    ALLOCATE(xx(ndim))
    xx = x0
    DO i = 1,ngrid(idim)
       pos(idim)=i
       IF(idim/=ndim-ndw+1) THEN
          IF (PRESENT(array)) THEN
             CALL fes_write(idim-1, fes,  pos,  ndim, ngrid, dp_grid,&
                  & x0, ndw, l_fes_int, array)
          ELSE
             CALL fes_write(idim-1, fes,  pos,  ndim, ngrid, dp_grid,&
                  & x0, ndw, l_fes_int)
          END IF
       ELSE
          IF (PRESENT(array)) THEN
             ind=1
             np=ngrid(ndim)*ngrid(ndim-1)*ngrid(ndim-2)
             DO is=1, ndw
                itt=1
                DO it=1,is-1
                   itt=itt*ngrid(ndim-it)
                END DO
                ind=ind+(pos(ndim-is+1)-1)*itt
             END DO
             IF (ind >np) STOP "something wrong in indexing .."
          END IF
          pnt=point_no_pbc(pos,ngrid,ndim)
          xx = x0 + dp_grid * ( pos - 1 )
          dimval = PRODUCT(ngrid(1:ndim-ndw))
          IF (l_fes_int) THEN
             IF (PRESENT(array)) THEN
                array(ind)=MINVAL(-fes(pnt:pnt+dimval-1))
             ELSE
                WRITE(123,'(10f20.10)')(xx(id),id=ndim,ndim-ndw+1,-1),MINVAL(-fes(pnt:pnt+dimval-1))
             END IF
          ELSE
             sum_fes =0.0_dp
             dvol = PRODUCT(dp_grid(1:ndim-ndw))
             DO is = pnt,pnt+dimval-1
                sum_fes = sum_fes + fes(is)*dvol
             END DO
             IF (PRESENT(array)) THEN
                array(ind)=-sum_fes
             ELSE
                WRITE(123,'(10f20.10)')(xx(id),id=ndim,ndim-ndw+1,-1), -sum_fes
             END IF
          END IF
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
  RECURSIVE SUBROUTINE fes_only_write(idim, fes,  pos,  ndim, ngrid, dp_grid, ndw, l_fes_int)
    INTEGER, INTENT(IN)                      :: idim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, DIMENSION(:), POINTER           :: pos
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid
    INTEGER, INTENT(IN)                      :: ndw
    LOGICAL, INTENT(IN)                      :: l_fes_int

    INTEGER                                  :: dimval, i, is, pnt
    REAL(KIND=dp)                            :: dvol, sum_fes

    DO i = 1,ngrid(idim)
       pos(idim)=i
       IF(idim/=ndim-ndw+1) THEN
          CALL fes_only_write(idim-1, fes,  pos,  ndim, ngrid, dp_grid, ndw, l_fes_int)
       ELSE        
          pnt=point_no_pbc(pos,ngrid,ndim)
          dimval = PRODUCT(ngrid(1:ndim-ndw))
          IF (l_fes_int) THEN
             WRITE(123,'(1f12.5)')MINVAL(-fes(pnt:pnt+dimval-1))
          ELSE
             sum_fes =0.0_dp
             dvol = PRODUCT(dp_grid(1:ndim-ndw))
             DO is = pnt,pnt+dimval-1
                sum_fes = sum_fes + fes(is)*dvol
             END DO
             WRITE(123,'(1f12.5)')-sum_fes
          END IF
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
  SUBROUTINE fes_min(fes, ndim, iperd, ngrid, dp_grid, x0, ndw)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: iperd, ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid, x0
    INTEGER, INTENT(IN)                      :: ndw

    INTEGER                                  :: i, id, iter, j, k, nacc, &
                                                ntrials, pnt
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: history
    INTEGER, DIMENSION(:), POINTER           :: pos, pos0
    INTEGER, DIMENSION(ndim)                 :: Dpos, ntrust
    LOGICAL                                  :: do_save
    REAL(KIND=dp)                            :: fes_now, fes_old, norm_dx, &
                                                resto
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dx, rnd, xx

    IF (ndw/=ndim) STOP "Not implemented for projected FES!"

    ntrust  = ngrid/10
    ntrials = PRODUCT (ngrid)
    WRITE(6,'(A,10I6)',ADVANCE="no")"FES| Trust hyper-radius ",ntrust
    WRITE(6,'(A,10F12.6)')" which is equivalent to: ",ntrust*dp_grid

    ALLOCATE(xx(ndim),dx(ndim),pos0(ndim),rnd(ndim),pos(ndim))
    ALLOCATE(history(ndim,ntrials))
    history = 0
    nacc    = 0
    Trials: DO j = 1, ntrials
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

       ! Loop over the frame points unless it is periodic
       DO k = 1, ndim
          IF ((iperd(k)==0).AND.(pos0(k)<ntrust(k)))          CYCLE Trials
          IF ((iperd(k)==0).AND.(pos0(k)>ngrid(k)-ntrust(k))) CYCLE Trials
       END DO

       ! Evaluate position and derivative
       pos = pos0
       xx  = x0 + dp_grid * ( pos - 1 )
       dx = derivative (fes, pos, iperd, ndim, ngrid, dp_grid)

       ! Integrate till derivative is small enough..
       pnt     = point_no_pbc(pos,ngrid,ndim)
       fes_now = -fes(pnt)
       fes_old = HUGE(0.0_dp)

       i = 1
       DO WHILE ((i<=100).OR.(fes_now < fes_old))
          fes_old = fes_now
          !WRITE(10+j,'(10f20.10)')(xx(id),id=ndim,1,-1),-fes(pnt)  

          norm_dx = SQRT(DOT_PRODUCT(dx,dx))
          IF (norm_dx == 0.0_dp) EXIT ! It is in a really flat region
          xx  = xx - MIN(0.1_dp,norm_dx)*dx/norm_dx
          ! Re-evaluating pos
          pos = CEILING((xx-x0)/dp_grid)+1
          CALL pbc(pos, iperd, ngrid, ndim)

          ! Incremental pos
          dx      = derivative (fes, pos, iperd, ndim, ngrid, dp_grid)
          pnt     = point_no_pbc(pos,ngrid,ndim)
          fes_now = -fes(pnt)
          i   = i + 1
       END DO
       iter = i

       ! Compare with the available minima and if they are the same skip
       ! saving this position..
       do_save = fes(pnt)>=1.0E-3_dp
       IF (do_save) THEN        
          DO i = 1, nacc
             Dpos = pos - history(:,i)
             norm_dx  = DOT_PRODUCT(Dpos,Dpos)
             IF ((SQRT(REAL(norm_dx,KIND=dp))<=MAXVAL(ntrust)).OR.(fes(pnt)<1.0E-3_dp)) THEN 
                do_save = .FALSE.
                EXIT
             END IF
          END DO
       END IF
       IF (do_save) THEN
          pnt = point_no_pbc(pos,ngrid,ndim)
          xx  = x0+dp_grid*(pos-1)
          WRITE(6,'(A,5F12.6)',ADVANCE="NO") "FES| Minimum found (",(xx(id),id=ndim,ndim-ndw+1,-1)
          WRITE(6,'(A,F12.6,A,I6)')" ). FES value = ",-fes(pnt)," Hartree. Number of Iter: ",iter
          nacc = nacc + 1
          history (:,nacc) = pos
       END IF
    END DO Trials
    WRITE(6,'(A,I6,A)') "FES| Number of Minimum found: ",nacc,"."

    DEALLOCATE(xx,dx,pos0,rnd,pos)
    DEALLOCATE(history)

  END SUBROUTINE fes_min

! *****************************************************************************
!> \brief Finds path between two points (a) and (b)
!> \par History
!>      12.2010 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fes_path(fes, ndim, ngrid, dp_grid, iperd, x0, ndw, mep_input_data, l_int)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid
    INTEGER, DIMENSION(:), POINTER           :: iperd
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    INTEGER, INTENT(IN)                      :: ndw
    TYPE(mep_input_data_type), INTENT(IN)    :: mep_input_data

    INTEGER                                  :: i, id, irep, iter, nf, &
                                                nreplica, ns, pnt
    INTEGER, DIMENSION(:), POINTER           :: ipos
    LOGICAL                                  :: converged, l_int
    REAL(KIND=dp)                            :: avg1, avg2, diff, ene, &
                                                norm_dx, xx0, yy0
    REAL(KIND=dp), DIMENSION(:), POINTER     :: davg1, davg2, dxx, dyy, &
                                                fes_rep, tang, xx, yy
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dx, pos, pos_old

    IF (ndw/=ndim) STOP "Not implemented for projected FES!"
    nreplica = mep_input_data%nreplica
    ALLOCATE(xx(ndim),dx(ndim,nreplica),pos_old(ndim,nreplica),pos(ndim,nreplica),&
             ipos(ndim),fes_rep(nreplica),dxx(ndim),dyy(ndim),yy(ndim),davg1(ndim),&
             tang(ndim),davg2(ndim))


    IF (l_int) THEN
        id = 0
      DO i=ndim,ndim-ndw+1,-1
         id = id +1
         pos(i,1)        = mep_input_data%minima(id,1)
         pos(i,nreplica) = mep_input_data%minima(id,2)
      END DO
    
      ! Interpolate nreplica-2 points
      xx=(pos(:,nreplica)-pos(:,1))/REAL(nreplica-1,KIND=dp)
      DO irep=2,nreplica-1
       pos(:,irep)=pos(:,1)+xx(:)*REAL(irep-1,KIND=dp)
      END DO

    ELSE
      pos=mep_input_data%minima
    ENDIF


    ! Compute value and derivative in all replicas
    DO irep=1,nreplica
       ipos= FLOOR((pos(:,irep)-x0)/dp_grid)+1
       pnt = point_no_pbc(ipos,ngrid,ndim)
       dx(:,irep) = derivative (fes, ipos, iperd, ndim, ngrid, dp_grid)
       fes_rep(irep) = -fes(pnt)
    END DO

    ! Implement a simple elastic band method (Hamiltonian): definitely not the best
    ! method, but for such a simple task it should be more than enough
    converged=.FALSE.
    pos_old=pos
    iter = 0
    DO WHILE ((.NOT.converged).AND.(iter<=mep_input_data%max_iter))
       iter=iter+1
       avg1=0.0_dp
       ! compute average length (distance 1)
       DO irep=2,nreplica
          xx=pos(:,irep)-pos(:,irep-1)
          avg1=avg1+SQRT(DOT_PRODUCT(xx,xx))
       END DO
       avg1=avg1/REAL(nreplica-1, KIND=dp)
       
       avg2=0.0_dp
       ! compute average length (distance 2)
       DO irep=3,nreplica
          xx=pos(:,irep)-pos(:,irep-2)
          avg2=avg2+SQRT(DOT_PRODUCT(xx,xx))
       END DO
       avg2=avg2/REAL(nreplica-2, KIND=dp)

       ! compute energy and derivatives
       dx = 0.0_dp
       ene= 0.0_dp
       ns=1
       nf=nreplica
       DO irep=1,nreplica
          ! compute energy and map point replica irep
          ipos= FLOOR((pos(:,irep)-x0)/dp_grid)+1
          pnt = point_no_pbc(ipos,ngrid,ndim)
          fes_rep(irep) = -fes(pnt)
          IF ((irep==1).OR.(irep==nreplica)) CYCLE

          ! -------------------------------------------------------------
          ! compute non-linear elastic terms : including only 2-d springs
          ! -------------------------------------------------------------
          davg2=0.0_dp
          IF (irep<nf-1) THEN
             xx  =pos(:,irep)-pos(:,irep+2)
             xx0 =SQRT(DOT_PRODUCT(xx,xx))
             dxx =1.0_dp/xx0*xx
             ene=ene+0.25_dp*mep_input_data%kb*(xx0-avg2)**2
             davg2=davg2+dxx
          END IF

          IF (irep>ns+1) THEN
             xx  =pos(:,irep)-pos(:,irep-2)
             yy0 =SQRT(DOT_PRODUCT(xx,xx))
             dyy =1.0_dp/yy0*xx
             davg2=davg2+dyy
          END IF
          davg2=davg2/REAL(nreplica-2, KIND=dp)
          
          IF (irep<nf-1) THEN
             dx(:,irep)=dx(:,irep)+0.5_dp*mep_input_data%kb*(xx0-avg2)*(dxx-davg2)
          END IF
          IF (irep>ns+1) THEN
             dx(:,irep)=dx(:,irep)+0.5_dp*mep_input_data%kb*(yy0-avg2)*(dyy-davg2)
          END IF

          ! -------------------------------------------------------------
          ! Evaluation of the elastic term
          ! -------------------------------------------------------------
          xx  =pos(:,irep)-pos(:,irep+1)
          yy0 =SQRT(DOT_PRODUCT(xx,xx))
          dyy =1.0_dp/yy0*xx

          xx  =pos(:,irep)-pos(:,irep-1)
          xx0 =SQRT(DOT_PRODUCT(xx,xx))
          dxx = 1.0_dp/xx0*xx
          davg1=(dxx+dyy)/REAL(nreplica-1, KIND=dp)

          ene=ene+0.5_dp*mep_input_data%kb*(xx0-avg1)**2
          dx(:,irep)=dx(:,irep)+mep_input_data%kb*(xx0-avg1)*(dxx-davg1)+&
                                mep_input_data%kb*(yy0-avg1)*(dyy-davg1)

          ! Evaluate the tangent
          xx=pos(:,irep+1)-pos(:,irep)
          xx=xx/SQRT(DOT_PRODUCT(xx,xx))
          yy=pos(:,irep)-pos(:,irep-1)
          yy=yy/SQRT(DOT_PRODUCT(yy,yy))
          tang=xx+yy
          tang=tang/SQRT(DOT_PRODUCT(tang,tang))

          xx= derivative (fes, ipos, iperd, ndim, ngrid, dp_grid)
          dx(:,irep) = DOT_PRODUCT(dx(:,irep),tang)*tang+&
                       xx - DOT_PRODUCT(xx,tang)*tang
       END DO
       dx(:,1)=0.0_dp
       dx(:,nreplica)=0.0_dp

       ! propagate the band with a SD step
       diff=0.0_dp
       DO irep=1,nreplica
          ene=ene+fes_rep(irep)
          IF ((irep==1).OR.(irep==nreplica)) CYCLE

          norm_dx = SQRT(DOT_PRODUCT(dx(:,irep),dx(:,irep)))
          IF (norm_dx/=0.0_dp) THEN
            pos(:,irep)  = pos(:,irep) - MIN(0.1_dp,norm_dx)*dx(:,irep)/norm_dx
          END IF
          xx=pos(:,irep)-pos_old(:,irep)
          diff=diff+DOT_PRODUCT(xx,xx)
       END DO
       IF (SQRT(diff)<=0.001_dp) THEN
          converged=.TRUE.
       END IF
       pos_old=pos
       WRITE(*,*)"Iteration nr.",iter, SQRT(diff)
    END DO

    WRITE(*,*)"MEP saved on <mep.data> file."
    OPEN(10,FILE="mep.data",STATUS='unknown',FORM='formatted')
    DO irep=1,nreplica
       ! compute energy and derivative for each single point of the replica
       ipos= FLOOR((pos(:,irep)-x0)/dp_grid)+1
       pnt = point_no_pbc(ipos,ngrid,ndim)
       fes_rep(irep) = -fes(pnt)
       WRITE(10,*)irep, pos(:,nreplica-irep+1), fes_rep(nreplica-irep+1)
    END DO
    CLOSE(10)

    DEALLOCATE(xx,dx,pos,fes_rep,ipos,pos_old,yy,dyy,dxx,davg1,tang,davg2)
  END SUBROUTINE fes_path

! *****************************************************************************
!> \brief Dump FES with a GAUSSIAN cube format - Useful for multidimensional FES
!> \par History
!>      12.2013 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fes_cube_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw, l_fes_int, file)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                      :: idim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fes
    INTEGER, DIMENSION(:), POINTER           :: pos
    INTEGER, INTENT(IN)                      :: ndim
    INTEGER, DIMENSION(:), POINTER           :: ngrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dp_grid, x0
    INTEGER, INTENT(IN)                      :: ndw
    LOGICAL, INTENT(IN)                      :: l_fes_int
    CHARACTER(LEN=80)                        :: file

    INTEGER                                  ::  ix, iy, iz, i, ip, natoms, np, &
         mesh(3), nspin, ind, id(3), iix, iiy, iiz, ii
    INTEGER, DIMENSION(:), POINTER           :: izat

    CHARACTER(LEN=5)                         :: label, labelp
    CHARACTER(LEN=120)                       :: line
    
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: rho, rhot
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: xat
    REAL(KIND=dp)                            :: cell(3,3), cm(3), rt(3),&
         delta(3), dr(3), residual
    
    CALL init_physcon
    IF (ndw .GT. 3) THEN
       WRITE(6,*) 
       WRITE(6,*) 'ERROR: GAUSSIAN format can only handle FES on 3 CV !'
       STOP
    END IF

    OPEN(10,file=file,status='old')
    CALL get_val_res(unit=10,section="&SUBSYS",subsection="&CELL")
    READ(10,*)label,cell(1,1),cell(2,1),cell(3,1) 
    READ(10,*)label,cell(1,2),cell(2,2),cell(3,2)
    READ(10,*)label,cell(1,3),cell(2,3),cell(3,3)
    rt(1)=-(cell(1,1)/2._dp)
    rt(2)=-(cell(2,2)/2._dp)
    rt(3)=-(cell(3,3)/2._dp)

    WRITE(6,*) 'Dumping GAUSSIAN CUBE format' 
    WRITE(6,*) 'Cell vectors'
    WRITE(6,*) 

    residual = 0.0d0
    DO ix=1,3
       DO iy=ix+1,3
          residual = residual + cell(ix,iy)**2
       END DO
    END DO

    IF (residual .gt. 1.0d-6) THEN
       WRITE(6,*) 
       WRITE(6,*) 'ERROR: this progam can only handle orthogonal cells'
       WRITE(6,*) ' with vectors pointing in the X, Y and Z directions'
       STOP
    END IF

    WRITE(6,*) 
    WRITE(6,*) 'Cube grid mesh: ',ngrid(1),'x',ngrid(2),'x',ngrid(3)
    WRITE(6,*) 'Origin in:',rt
    WRITE(6,*) 
    
    DO ix=1,3
       dr(ix)=cell(ix,ix)/REAL(ngrid(ix),KIND=dp)
    END DO
    
    np = ngrid(1) * ngrid(2) * ngrid(3)    
    ALLOCATE(rho(np), rhot(np))
    CALL fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw, l_fes_int, rho)
    WRITE(6,*) 'Internal FES transfer completed!'

    ! translate cell
    DO ix=1,3
       delta(ix) = rt(ix)/dr(ix)
       id(ix) = delta(ix)
       delta(ix) = rt(ix) - id(ix) * dr(ix)
    END DO

    DO iz=1,ngrid(3)
       DO iy=1,ngrid(2)
          DO ix=1,ngrid(1)
             iix=ix+id(1)
             iiy=iy+id(2)
             iiz=iz+id(3)
             IF (iix .LT. 1) iix=iix+ngrid(1)
             IF (iiy .LT. 1) iiy=iiy+ngrid(2)
             IF (iiz .LT. 1) iiz=iiz+ngrid(3)
             IF (iix .GT. ngrid(1)) iix=iix-ngrid(1)
             IF (iiy .GT. ngrid(2)) iiy=iiy-ngrid(2)
             IF (iiz .GT. ngrid(3)) iiz=iiz-ngrid(3)
             
             IF (iix .LT. 1) STOP 'ix < 0'
             IF (iiy .LT. 1) STOP 'iy < 0'
             IF (iiz .LT. 1) STOP 'iz < 0'
             IF (iix .GT. ngrid(1)) STOP 'ix > cell'
             IF (iiy .GT. ngrid(2)) STOP 'iy > cell'
             IF (iiz .GT. ngrid(3)) STOP 'iz > cell'
             i=ix+(iy-1)*ngrid(1)+(iz-1)*ngrid(1)*ngrid(2)
             ii=iix+(iiy-1)*ngrid(1)+(iiz-1)*ngrid(1)*ngrid(2)
             rhot(ii)=rho(i)
          END DO
       END DO
    END DO
    
    REWIND(10)
    CALL get_val_res(unit=10,section="&SUBSYS",subsection="&COORD")
    natoms=0
    ALLOCATE(xat(1000,3))
    ALLOCATE(izat(1000))
    DO WHILE (.TRUE.)
       READ(10,'(A)') line
       IF (INDEX(line,'&END') /= 0) EXIT
       natoms=natoms+1
       READ(line,*)label,(xat(natoms,ix),ix=1,3)
       IF (natoms == SIZE(xat,1)) THEN
          CALL reallocate(xat,1,SIZE(xat,1)*2,1,3)
          CALL reallocate(izat,1,SIZE(izat)*2)
       END IF
       CALL uppercase(label)
       DO i=1, nelem
          labelp=ptable(i)%symbol
          CALL uppercase(labelp)
          IF (TRIM(label)==TRIM(labelp)) EXIT
       END DO
       IF (i == nelem+1) THEN
          WRITE(6,*)TRIM(label), "In line: ", line
          STOP "Element not recognized!"
       END IF
       izat(natoms)=i
    END DO
    CALL reallocate(xat,1,natoms,1,3)
    CALL reallocate(izat,1,natoms)
    
    DO i=1,natoms
       DO ix=1,3
          xat(i,ix)=xat(i,ix)+rt(ix)-delta(ix)
          IF (xat(i,ix) .LT. 0.0)         xat(i,ix)=xat(i,ix)+cell(ix,ix)
          IF (xat(i,ix) .GT. cell(ix,ix)) xat(i,ix)=xat(i,ix)-cell(ix,ix)
       END DO
    END DO
    
    WRITE(123,*) "FES on CUBE"
    WRITE(123,*) "created by fes in CP2K"
    WRITE(123,'(i5,4f12.6)') natoms, 0.0,0.0,0.0

    DO ix=1,3
       ii = ngrid(ix)
       WRITE(123,'(i5,4f12.6)') ii,(cell(ix,iy)/ii*bohr,iy=1,3)
    END DO

    DO i=1,natoms
       WRITE(123,'(i5,4f12.6)') izat(i),0.0,(xat(i,ix)*bohr,ix=1,3)
    END DO

    DO ix=1,ngrid(1)
       DO iy=1,ngrid(2)

        WRITE(123,'(6e13.5)')(rhot(ix+(iy-1)*ngrid(1)+(iz-1)*ngrid(1)&
             &*ngrid(2)),iz=1,ngrid(3)) 
     END DO
  END DO
  DEALLOCATE(xat, rho, rhot)

  END SUBROUTINE fes_cube_write

END MODULE graph_methods
