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
PROGRAM fes_compute

  USE kinds,                           ONLY: dp,&
                                             sp,&
                                             default_string_length
  USE machine,                         ONLY: m_getlog,&
                                             m_hostnm,&
                                             m_getarg,&
                                             m_iargc
  USE mathconstants,                   ONLY: pi
  
  IMPLICIT NONE
  
  INTERFACE
! *****************************************************************************
     SUBROUTINE get_val_res(unit, section, keyword, subsection, i_val)
       IMPLICIT NONE
       INTEGER, INTENT(IN)            :: unit
       CHARACTER (len=*)              :: section 
       CHARACTER (len=*),OPTIONAL     :: subsection 
       CHARACTER (len=*),OPTIONAL     :: keyword
       INTEGER, OPTIONAL, INTENT(OUT) :: i_val
     END SUBROUTINE get_val_res

! *****************************************************************************
     RECURSIVE SUBROUTINE fes_compute_low(idim, nn, fes, gauss, ind, ind0, nfes, ndim, ngauss, ngrid, iperd)
       USE kinds,                           ONLY: dp,&
                                                  sp       
       IMPLICIT NONE
       INTEGER, INTENT(in) :: idim, nfes, ngauss, ndim
       INTEGER, DIMENSION(:) :: nn, ind, ind0, iperd
       INTEGER, DIMENSION(:), POINTER :: ngrid
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes
       REAL(KIND=dp), DIMENSION(:,:), POINTER :: gauss
     END SUBROUTINE fes_compute_low

! *****************************************************************************
     RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw)
       USE kinds,                           ONLY: dp,&
                                                  sp
       IMPLICIT NONE
       INTEGER, INTENT(in) :: idim, ndim, ndw
       INTEGER, DIMENSION(:), POINTER :: pos
       INTEGER, DIMENSION(:), POINTER :: ngrid
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0
     END SUBROUTINE fes_write

! *****************************************************************************
     RECURSIVE SUBROUTINE fes_only_write(idim, fes,  pos,  ndim, ngrid, ndw)
       USE kinds,                           ONLY: dp,&
                                                  sp       
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ndw
       INTEGER, INTENT(in) :: idim, ndim
       INTEGER, DIMENSION(:), POINTER :: ngrid
       INTEGER, DIMENSION(:), POINTER :: pos
       
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes
     END SUBROUTINE fes_only_write

! *****************************************************************************
     SUBROUTINE fes_min(fes,  ndim, ngrid, dp_grid, x0, ndw)
       USE kinds,                           ONLY: dp,&
                                                  sp
       IMPLICIT NONE
       INTEGER, INTENT(in) :: ndim, ndw
       INTEGER, DIMENSION(:), POINTER :: ngrid
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0
     END SUBROUTINE fes_min

! *****************************************************************************
     SUBROUTINE fes_path(fes,  ndim, ngrid, dp_grid, x0, ndw, minima)
       USE kinds,                           ONLY: dp,&
                                                  sp
       IMPLICIT NONE
       REAL(KIND=dp), DIMENSION(:,:), POINTER :: minima
       INTEGER, INTENT(in) :: ndim, ndw
       INTEGER, DIMENSION(:), POINTER :: ngrid
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0
     END SUBROUTINE fes_path

  END INTERFACE

  INTEGER :: i, ix ,id, it, ip, stat, itmp, iw
  INTEGER :: ndim, ndw, nt_p, nt, nh, nf, nwr, ncount, nprd
  INTEGER :: argcount, ngauss, nfes, coor
  INTEGER, POINTER :: idw(:), iperd(:), iprd(:)
  INTEGER, POINTER :: i_map(:)
  INTEGER, POINTER :: ngrid(:), tmp(:)
  INTEGER, POINTER :: ind(:), inds(:), nn(:,:), nn_max(:)

  REAL(KIND=dp) :: dum, ss, diff
  REAL(KIND=dp) :: x0w(3),xfw(3)
  REAL(KIND=dp), POINTER ::dp_cut(:)
  REAL(KIND=dp), POINTER :: ss0(:,:)
  REAL(KIND=dp), POINTER :: delta_s(:,:)
  REAL(KIND=dp), POINTER :: gauss(:,:),fes(:)
  REAL(KIND=dp), POINTER :: ww(:), minima(:,:)
  REAL(KIND=dp), POINTER :: dp_grid(:),x0(:),xf(:), tmpr(:)
  REAL(KIND=dp) :: dp2,  eps_cut

  CHARACTER(LEN=80) :: wq_char, file, out

  LOGICAL :: fix, lstride
  LOGICAL :: l_grid, l_dp, l_orac, l_cp2k, l_math, l_fmin, l_pmin
  CHARACTER(LEN=default_string_length) :: per_label
  INTRINSIC :: random_seed, random_number

  ! Initialize variables
  nprd    = 0
  ndim    = 1
  ndw     = 1
  nt_p    = 9999999
  eps_cut = 1e-6
  file    = 'HILLS'
  out     = 'fes.dat'
  fix     = .FALSE.
  lstride = .FALSE.
  l_grid  = .FALSE.
  l_dp    = .FALSE.
  l_orac  = .FALSE.
  l_cp2k  = .FALSE.
  l_math  = .FALSE.
  l_fmin  = .FALSE.
  l_pmin  = .FALSE.
  iw      = 6

  argcount = m_iargc()
  IF(argcount==0)THEN
     WRITE(iw,*)'USAGE:'
     WRITE(iw,*)'graf  '
     WRITE(iw,*)'[-ngrid  50 .. ..]   (Mesh dimension. Default :: 100)'
     WRITE(iw,*)'[-dp   0.05 .. ..]   (Alternative to -ngrid, allows the specification of the mesh dx)'
     WRITE(iw,*)'[-ndim  3        ]   (Number of collective variables NCV)'
     WRITE(iw,*)'[-ndw  1 3  ..   ]   (CVs for the free energy surface)'
     WRITE(iw,*)'[-periodic 2 3 ..]   (CVs with periodic boundary conditions (-pi,pi] )'
     WRITE(iw,*)'[-stride 10      ]   (How often the FES is written)'
     WRITE(iw,*)'[-fix   1.1 .. ..]   (Define the region for the FES)'
     WRITE(iw,*)'                     (If omitted this is automatically calculated)'
     WRITE(iw,*)'[-cutoff 2.      ]   (The hills are cutoffed at 2)'
     WRITE(iw,*)'[-file   filename]'
     WRITE(iw,*)'[-out    filename]'
     WRITE(iw,*)'[-orac]              (If energies are written in orac intern units)'
     WRITE(iw,*)'[-cp2k]              (Specify if a CP2K restart file is provided)'
     WRITE(iw,*)'[-mathlab]           (File storing FES in Mathlab format. Default format Gnuplot)'
     WRITE(iw,*)'[-find-minima]       (Tries to finds all minima in the computed FES)'
     WRITE(iw,*)'[-find-path]         (Finds MEP between all minima (found) in the computed FES)'
     WRITE(iw,*)'[-point-a]           (Specifies point (a) when using -find-path option)'
     WRITE(iw,*)'[-point-b]           (Specifies point (b) when using -find-path option)'
     WRITE(iw,*)''
     WRITE(iw,*)'DEFAULT OUTPUT: fes.dat'
     WRITE(iw,*)''
     STOP  "Please provide arguments to run FES!"
  ENDIF
  
  DO i=1,argcount
     CALL M_GETARG(i,wq_char)
 
     IF (INDEX(wq_char,'-file').NE.0)THEN
        CALL M_GETARG(i+1,wq_char)
        READ(wq_char,*)file
     ENDIF

     IF (INDEX(wq_char,'-out').NE.0)THEN
        CALL M_GETARG(i+1,wq_char)
        READ(wq_char,*)out
     ENDIF

     IF (INDEX(wq_char,'-ndim').NE.0)THEN
        CALL M_GETARG(i+1,wq_char)
        READ(wq_char,*)ndim
     ENDIF

     IF (INDEX(wq_char,'-stride').NE.0)THEN
        CALL M_GETARG(i+1,wq_char)
        READ(wq_char,*)nt_p
        lstride=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-cutoff').NE.0)THEN
        CALL M_GETARG(i+1,wq_char)
        READ(wq_char,*)eps_cut
     ENDIF

     IF (INDEX(wq_char,'-orac').NE.0)THEN
        l_orac=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-cp2k').NE.0)THEN
        l_cp2k=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-find-minima').NE.0)THEN
        l_fmin=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-find-path').NE.0)THEN
        l_pmin=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-mathlab').NE.0)THEN
        l_math=.TRUE.
     ENDIF
  END DO
  IF (l_cp2k) l_orac = .FALSE.

  ! Initializing random numbers
  CALL RANDOM_SEED()
  CALL RANDOM_NUMBER(dum)

  ! Basic Allocation
  ndw = ndim
  ALLOCATE(ngrid(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE(dp_grid(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE(idw(ndw),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE(iperd(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE(iprd(nprd),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  DO i = 1, ndim
     idw(i)   = i
     iperd(i) = 0
  END DO

  DO i=1,argcount
     CALL M_GETARG(i,wq_char)

     IF (INDEX(wq_char,'-ndw').NE.0)THEN
        DEALLOCATE(idw,stat=stat)
        IF (stat/=0) STOP "Allocation Error"

        ndw=0
        ndw_loop: DO ix=i+1,argcount
           CALL M_GETARG(ix,wq_char)
           IF(INDEX(wq_char,'-').EQ.0)THEN
              ndw=ndw+1
           ELSE
              EXIT ndw_loop
           ENDIF
        ENDDO ndw_loop

        ALLOCATE(idw(ndw),stat=stat)
        IF (stat/=0) STOP "Allocation Error"

        DO id=1,ndw
           CALL M_GETARG(i+id,wq_char)
           READ(wq_char,*)idw(id)
        ENDDO
     ENDIF

     IF (INDEX(wq_char,'-periodic').NE.0)THEN
        nprd=0
        nprd_loop: DO ix=i+1,argcount
           CALL M_GETARG(ix,wq_char)
           IF(INDEX(wq_char,'-').EQ.0)THEN
              nprd=nprd+1
           ELSE
              EXIT nprd_loop
           ENDIF
        ENDDO nprd_loop

        DEALLOCATE(iprd,stat=stat)
        IF (stat/=0) STOP "Allocation Error"
        ALLOCATE(iprd(nprd),stat=stat)
        IF (stat/=0) STOP "Allocation Error"

        DO id=1,nprd
           CALL M_GETARG(i+id,wq_char)
           READ(wq_char,*)iprd(id)
        ENDDO
     ENDIF

     IF (INDEX(wq_char,'-ngrid').NE.0)THEN
        DO ix=1,ndim
           CALL M_GETARG(i+ix,wq_char)
           READ(wq_char,*)ngrid(ix)
           l_grid=.TRUE.
        END DO
     ENDIF

     IF (INDEX(wq_char,'-dp').NE.0)THEN
        DO ix=1,ndim
           CALL M_GETARG(i+ix,wq_char)
           READ(wq_char,*)dp_grid(ix)
           l_dp=.TRUE.
           l_grid=.FALSE.
        END DO
     END IF
     
     IF (INDEX(wq_char,'-fix').NE.0)THEN
        fix=.TRUE.
        DO id=1,ndw
           CALL M_GETARG(i+2*(id-1)+1,wq_char)
           READ(wq_char,*)x0w(id)
           CALL M_GETARG(i+2*(id-1)+2,wq_char)
           READ(wq_char,*)xfw(id)
        ENDDO
     ENDIF
  ENDDO

  IF (l_pmin) THEN
     ALLOCATE(minima(ndw,2))
     minima = HUGE(0.0_dp)
     ! Read for starting point (a) and (b)
     DO i=1,argcount
        CALL M_GETARG(i,wq_char)
        
        IF (INDEX(wq_char,'-point-a').NE.0)THEN
           DO id=1,ndw
              CALL M_GETARG(i+id,wq_char)
              READ(wq_char,*)minima(id,1)
           ENDDO
        ENDIF

        IF (INDEX(wq_char,'-point-b').NE.0)THEN
           DO id=1,ndw
              CALL M_GETARG(i+id,wq_char)
              READ(wq_char,*)minima(id,2)
           ENDDO
        ENDIF
     END DO
     IF (ANY(minima==HUGE(0.0_dp))) STOP "-find-path requires the specification of -point-a and -point-b !"
  ELSE
     ALLOCATE(minima(0,0))
  END IF

  !  Defines the order of the collectiv var.: first the "wanted" ones, then the others
  ALLOCATE(i_map(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error" 
  i_map = 0

  DO id=1,ndw
     i_map(idw(id))=id
  ENDDO
  ix=ndw
  DO id=1,ndim
     IF(i_map(id)==0)THEN
        ix=ix+1
        i_map(id)=ix
     ENDIF
  ENDDO
  i_map=ndim-i_map+1

  ! Tag the periodic COLVAR according the new internal order
  DO id=1,nprd
     DO i = 1, ndim
        IF (i_map(i)==iprd(id)) EXIT
     END DO
     IF (i==ndim+1) STOP "Error assigning the periodicity of COLVAR" 
     iperd(i)=1
  END DO

  ! Grid size
  IF(l_grid) THEN
     ALLOCATE(tmp(ndim),stat=stat)
     IF (stat/=0) STOP "Allocation Error"
     tmp=ngrid
     DO i=1,ndim
        ngrid(i_map(i))=tmp(i)
     END DO
     DEALLOCATE(tmp,stat=stat)
     IF (stat/=0) STOP "Allocation Error" 
  ELSE
     ngrid=100
  END IF

  WRITE(iw,'(/,70("*"))')
  WRITE(iw,'("FES|",T7,A,/)')"Parsing file:   <"//TRIM(file)//">"

  OPEN(10,file=file,status='old')
  IF (l_cp2k) THEN
     CALL get_val_res(unit=10,section="&METADYN",keyword="NHILLS_START_VAL", i_val=nt)
  ELSE
     nt=0
     DO WHILE (.TRUE.)
        READ(10,*,END=100,ERR=100)dum
        nt=nt+1
     END DO
100  REWIND(10)
  END IF
  
  ALLOCATE( x0(ndim)         , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( xf(ndim)         , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( ss0(ndim,nt)     , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( delta_s(ndim,nt) , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( ww(nt)           , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( ind(ndim)        , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( inds(ndim)       , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( nn(ndim,nt)      , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( nn_max(ndim)     , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( dp_cut(ndim)     , stat=stat )
  IF (stat/=0) STOP "Allocation Error"

  IF (l_cp2k) THEN
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_POS")
     DO i = 1, nt
        READ(10,*)(ss0(i_map(id),i),id=1,ndim)
     END DO
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_AMPL")
     DO i = 1, nt
        READ(10,*)(delta_s(i_map(id),i),id=1,ndim)
     END DO
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_HEIGHT")
     DO i = 1, nt
        READ(10,*)ww(i)
     END DO
  ELSE
     DO i = 1, nt
        READ(10,*)dum,(ss0(i_map(id),i),id=1,ndim),(delta_s(i_map(id),i),id=1,ndim),ww(i)
     END DO
  END IF
  CLOSE(10)

  ! ORAC conversion factor
  IF(l_orac) ww = ww * 10000._dp / 4.187_dp

  ! Setting up the limit of definitions for the several colvars
  DO id=1,ndim
     x0(id) = HUGE(1.0_dp)
     xf(id) =-HUGE(1.0_dp)
  ENDDO
  IF(fix) THEN
     DO it=1,nt
        DO id=1,ndim-ndw
           x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
           xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
        ENDDO
     ENDDO
     it=0
     DO id=ndim,ndim-ndw+1,-1
        it=it+1
        x0(id)=x0w(it)
        xf(id)=xfw(it)
     ENDDO
  ELSE
     DO it=1,nt
        DO id=ndim,1,-1
           x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
           xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
        ENDDO
     ENDDO
  ENDIF

  DO id=ndim,1,-1
     IF (iperd(id)==1) THEN
        x0(id)=-pi
        xf(id)= pi
     END IF
  END DO

  IF(l_dp)THEN
     ALLOCATE(tmpr(ndim))
     tmpr=dp_grid
     DO i=1,ndim
        dp_grid(i_map(i))=tmpr(i)
     END DO
     DEALLOCATE(tmpr)
     ngrid=INT((xf-x0)/dp_grid)+1
  ELSE
     dp_grid=(xf-x0)/DBLE(ngrid-1)
  END IF
  
  WRITE(iw,'(70("*"))')
  WRITE(iw,'("FES|",T7,A,/)')"Parameters for FES:"
  WRITE(iw,'("FES|",T7,A15,5x,i7)')"NDIM         ::",ndim
  WRITE(iw,'("FES|",T7,A15,5x,i7)')"NWD          ::",ndw
  WRITE(iw,'("FES|",T7,A15,5x,i7)')"HILLS        ::",nt
  it=ndim
  DO i=1,ndim
     per_label = ""
     IF (iperd(i)/=0) per_label = "(PERIODIC COLVAR)"
     WRITE(iw,'("FES|",T7,"COLVAR # ",i3," ::",5x,"(",f7.3," ,",f7.3,")",5X,A)')i_map(i),x0(i),xf(i),TRIM(per_label)
  END DO
  WRITE(iw,'("FES|",T7,a15,5x,7i7)'   )"NGRID        ::",(ngrid(id),id=ndim,ndim-ndw+1,-1)
  WRITE(iw,'("FES|",T7,a15,5x,5f7.3)' )"DX           ::",(dp_grid(id),id=ndim,ndim-ndw+1,-1)
  WRITE(iw,'("FES|",T7,a15,5x,g10.5)' )"CUTOFF       ::",eps_cut
  WRITE(iw,'(70("*"),/)')

  nn_max = 0
  DO i = 1, nt
     dp_cut  = SQRT(LOG(ABS(ww(i))/eps_cut))*2.0_dp*delta_s(:,i) 
     nn(:,i) = INT(dp_cut/dp_grid)
     ww(i)   = ww(i)**(1.0_dp/DBLE(ndim))
  END DO

  nn_max = MAXVAL(nn,DIM=2)
  ngauss = MAXVAL(nn_max) * 2 + 1
  nfes   = PRODUCT(ngrid)

  ALLOCATE(gauss(-MAXVAL(nn_max):MAXVAL(nn_max),ndim))
  ALLOCATE(fes(nfes))
  fes=0.0_dp

  nh=1
  nf=MIN(nh+nt_p-1,nt)
  
  IF (lstride)THEN 
     nwr=nt_p
  ELSE
     nwr=INT(nt/10)+1
  END IF

  ncount = 0
  WRITE(iw,'(/,"FES|",T7,A)') "Computing Free Energy Surface"
  Stride : DO WHILE (nh <= nt)
     Hills : DO it=nh,nf
        ind=INT((ss0(:,it)-x0)/dp_grid) + 1
        gauss=0.0_dp
        
        DO i=1,ndim
           coor = ind(i) - nn(i,it) - 1
           ss  = x0(i) + coor * dp_grid(i) - dp_grid(i)
           DO ip=-nn(i,it),nn(i,it)
              coor = coor + 1
              ss = ss + dp_grid(i)
              IF (iperd(i)==0) THEN
                 IF (coor .GT. ngrid(i)) CYCLE
                 IF (coor .LT. 1) CYCLE
              END IF
              diff = ss-ss0(i,it)
              dp2=(diff/delta_s(i,it))**2
              gauss(ip,i)=ww(it)*EXP(-0.5_dp*dp2)
           END DO
        END DO
        inds = ind
        CALL fes_compute_low(ndim,nn(:,it),fes,gauss,ind,inds,nfes,ndim,ngauss,ngrid,iperd)
        
        IF(.NOT. lstride .AND. MOD(it,nwr)==0)THEN
           WRITE(iw,'("FES|",T7,a,i4,a2)') "Mapping Gaussians ::",INT(10*ANINT(10.*it/nt))," %"
        ELSEIF(.NOT. lstride .AND. it==nt)THEN
           WRITE(iw,'("FES|",T7,a,i4,a2)') "Mapping Gaussians ::",INT(10*ANINT(10.*it/nt))," %"
        END IF
     END DO Hills

     IF (lstride) THEN
        ncount = ncount+1
        WRITE(iw,'("FES|",T7,a13,i5," |-| Gaussians from ",i6," to",i6)') "Done frame ::",ncount,nh,nf
        IF(ncount<10) THEN
           WRITE(file,'("fes.dat.",i1)')ncount
        ELSEIF(ncount<100) THEN
           WRITE(file,'("fes.dat.",i2)')ncount
        ELSE
           WRITE(file,'("fes.dat.",i3)')ncount
        END IF
        OPEN(123,file=file)
        ind   = 1
        CALL fes_only_write(ndim, fes,  ind,  ndim, ngrid, ndw)
        CLOSE(123)
     END IF

     nh=nh+nt_p
     nf=MIN(nh+nt_p-1,nt)
  END DO Stride
  DEALLOCATE(gauss)

  WRITE(iw,'("FES|",T7,A)') "Dumping FES structure in file: < "//TRIM(out)//" >"
  OPEN(123,file=out)
  ix=0
  IF (l_math) WRITE(123,'(10g12.5)')(ngrid(id),id=ndim,ndim-ndw+1,-1),ix
  ind   = 1   
  CALL fes_write(ndim, fes, ind, ndim, ngrid, dp_grid, x0, ndw)
  CLOSE(123)
 
  ! If requested find minima
  IF (l_fmin) CALL fes_min(fes, ndim, ngrid, dp_grid, x0, ndw)

  ! If requested find path
  IF (l_pmin) CALL fes_path(fes, ndim, ngrid, dp_grid, x0, ndw, minima)
     
  ! Free memory
  DEALLOCATE(ngrid,stat=stat)
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(dp_grid,stat=stat)
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(idw,stat=stat)
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(iperd,stat=stat)
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(x0,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(xf,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(ss0,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(delta_s,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(ww,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(ind,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(inds,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(nn,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(nn_max,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(dp_cut,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(i_map,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(fes,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(iprd,stat=stat )
  IF (stat/=0) STOP "Deallocation Error"
  DEALLOCATE(minima,stat=stat)
  IF (stat/=0) STOP "Deallocation Error"

  ! Terminate FES
  WRITE(iw,'(/,A,/)') "FES| NORMAL FES TERMINATION."
END PROGRAM fes_compute

! *****************************************************************************
!> \brief Efficiently map the gaussians on the grid
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
RECURSIVE SUBROUTINE fes_compute_low(idim, nn, fes, gauss, ind, ind0, nfes, ndim, ngauss, ngrid, iperd)
  USE kinds,                           ONLY: dp,&
                                             sp       
  IMPLICIT NONE
  INTEGER :: i,pnt,j,k
  INTEGER, INTENT(in) :: idim, nfes, ngauss, ndim
  INTEGER, DIMENSION(:) :: nn, ind, ind0, iperd
  INTEGER, DIMENSION(:), POINTER :: pos , ll, ngrid
  
  REAL(KIND=dp), DIMENSION(:), POINTER   :: fes
  REAL(KIND=dp), DIMENSION(:,:), POINTER :: gauss
  REAL(KIND=dp) :: prod

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
        pnt  = point(pos,iperd)
        prod = 1.0_dp
        DO j=1,ndim
           ll(j) = pos(j) - ind0(j)
           prod  = prod*gauss(ll(j),j)
        END DO
        fes(pnt)=fes(pnt)+prod
     END IF
  END DO
  DEALLOCATE(pos,ll)
  
CONTAINS
  
  INTEGER FUNCTION point (pos,iperd) RESULT(pnt)
    INTEGER, DIMENSION (:) :: pos, iperd
    INTEGER :: idim, lpnt
    
    idim=1
    pnt=pos(idim)
    IF (iperd(idim) == 1) THEN
       lpnt = pos(idim)
       lpnt = 1000*ngrid(idim)+lpnt
       lpnt = MOD(lpnt,ngrid(idim))
       IF (lpnt==0) lpnt=ngrid(idim)
       pnt = lpnt
    END IF
    DO idim=2,ndim
       lpnt = pos(idim)
       IF (iperd(idim) == 1) THEN
          lpnt = 1000*ngrid(idim)+lpnt
          lpnt = MOD(lpnt,ngrid(idim))
          IF (lpnt==0) lpnt=ngrid(idim)
       END IF
       pnt=pnt+(lpnt-1) * PRODUCT(ngrid(1:idim-1))
    END DO
  END FUNCTION point
  
END SUBROUTINE fes_compute_low

! *****************************************************************************
!> \brief Parser informations from the cp2k input/restart
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
SUBROUTINE get_val_res(unit, section, keyword, subsection, i_val)
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: unit
  CHARACTER (len=*)              :: section 
  CHARACTER (len=*),OPTIONAL     :: subsection 
  CHARACTER (len=*),OPTIONAL     :: keyword
  INTEGER, OPTIONAL, INTENT(OUT) :: i_val
  ! Variables
  INTEGER :: stat, my_ind
  CHARACTER (len=512)        :: line

  REWIND(unit)
  CALL search(unit,TRIM(section),line,stat=stat)
  
  IF (stat/=0) THEN 
     WRITE(6,*)"Pattern: "//TRIM(section)//" not found in input file!"
     STOP "Search failed!"
  END IF

  IF (PRESENT(keyword)) THEN
     CALL search(unit,TRIM(keyword),line,stat)
     my_ind=INDEX(line,TRIM(keyword))+LEN_TRIM(keyword)+1
     IF (PRESENT(i_val)) READ(line(my_ind:),*)i_val
  END IF
     
  IF (PRESENT(subsection)) THEN
     CALL search(unit,TRIM(subsection),line,stat)
  END IF

CONTAINS
  
! *****************************************************************************
  SUBROUTINE search(unit, key, line, stat)
    IMPLICIT NONE
    INTEGER, INTENT(in)            :: unit
    CHARACTER(LEN=*), INTENT(IN)   :: key
    CHARACTER(LEN=512),INTENT(OUT) :: line
    INTEGER, INTENT(out)           :: stat
    
    stat = 99
    DO WHILE (.TRUE.)
       READ(unit,'(A)',ERR=100,END=100) line
       IF (INDEX(line,TRIM(key)) /= 0) THEN 
          stat = 0
          EXIT 
       END IF
    END DO
100 CONTINUE
  END SUBROUTINE search

END SUBROUTINE get_val_res

! *****************************************************************************
!> \brief Writes the FES on the file
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw)
  USE kinds,                           ONLY: dp,&
                                             sp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndw
  INTEGER :: i,pnt,id,dimval
  INTEGER, INTENT(in) :: idim, ndim
  INTEGER, DIMENSION(:), POINTER :: pos, ngrid
  
  REAL(KIND=dp), DIMENSION(:), POINTER   :: xx
  REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0

  ALLOCATE(xx(ndim))
  xx = x0
  DO i = 1,ngrid(idim)
     pos(idim)=i
     IF(idim/=ndim-ndw+1) THEN
        CALL fes_write(idim-1, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw)
     ELSE        
        pnt=point(pos)
        xx = x0 + dp_grid * ( pos - 1 )
        dimval = PRODUCT(ngrid(1:ndim-ndw)) 
        WRITE(123,'(10f20.10)')(xx(id),id=ndim,ndim-ndw+1,-1),MINVAL(-fes(pnt:pnt+dimval-1))
     END IF     
  END DO
  DEALLOCATE(xx)
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE fes_write

! *****************************************************************************
!> \brief Writes the FES on the file when stride is requested
!> \par History
!>      03.2006 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
RECURSIVE SUBROUTINE fes_only_write(idim, fes,  pos,  ndim, ngrid, ndw)
  USE kinds,                           ONLY: dp,&
                                             sp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndw
  INTEGER :: i,pnt,dimval
  INTEGER, INTENT(in) :: idim, ndim
  INTEGER, DIMENSION(:), POINTER :: pos, ngrid
  
  REAL(KIND=dp), DIMENSION(:), POINTER   :: fes

  DO i = 1,ngrid(idim)
     pos(idim)=i
     IF(idim/=ndim-ndw+1) THEN
        CALL fes_only_write(idim-1, fes,  pos,  ndim, ngrid, ndw)
     ELSE        
        pnt=point(pos)
        dimval = PRODUCT(ngrid(1:ndim-ndw))
        WRITE(123,'(1f12.5)')MINVAL(-fes(pnt:pnt+dimval-1))
     END IF
  END DO

CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE fes_only_write

! *****************************************************************************
!> \brief Finds minima of the FES
!> \par History
!>      06.2009 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
SUBROUTINE fes_min(fes, ndim, ngrid, dp_grid, x0, ndw)
  USE kinds,                           ONLY: dp,&
                                             sp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndw
  INTEGER :: i,j,k,pnt,id,ntrials,nacc
  INTEGER, INTENT(in) :: ndim
  INTEGER, DIMENSION(:), POINTER :: pos, ngrid, pos0

  REAL(KIND=dp), DIMENSION(:), POINTER   :: xx, dx, rnd
  REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0
  REAL(KIND=dp)                          :: norm_dx, resto
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: history
  INTEGER, DIMENSION(ndim)               :: Dpos, ntrust
  LOGICAL                                :: do_save

  ALLOCATE(xx(ndim),dx(ndim),pos0(ndim),rnd(ndim),pos(ndim))
  ALLOCATE(history(ndim,ntrials))

  ntrust  = ngrid/10
  ntrials = PRODUCT (ngrid)
  WRITE(6,'(A,10I6)',ADVANCE="no")"FES| Trust hyper-radius ",ntrust
  WRITE(6,'(A,10F12.6)')" which is equivalent to: ",ntrust*dp_grid

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
        pnt = point(pos,ndim,ngrid)
        !WRITE(10+j,'(10f20.10)')(xx(id),id=ndim,1,-1),-fes(pnt)  
        
        norm_dx = SQRT(DOT_PRODUCT(dx,dx))
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
        pnt = point(pos,ndim,ngrid)
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

CONTAINS
  
  FUNCTION point (pos,ndim,ngrid) RESULT(pnt)
    IMPLICIT NONE
    INTEGER                           :: pnt
    INTEGER, DIMENSION(:), INTENT(IN) :: pos
    INTEGER, INTENT(IN)               :: ndim
    INTEGER, DIMENSION(:), INTENT(IN) :: ngrid
    INTEGER :: i

    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
  END FUNCTION point

  FUNCTION derivative(fes,pos0, ndim, ngrid, dp_grid) RESULT(der)
    USE kinds,                           ONLY: dp,&
                                               sp
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN)       :: pos0
    INTEGER, INTENT(IN)                     :: ndim
    INTEGER, DIMENSION(:), INTENT(IN)       :: ngrid
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: dp_grid, fes
    REAL(KIND=dp), DIMENSION(ndim)          :: der
    INTEGER, DIMENSION(:), ALLOCATABLE      :: pos
    INTEGER :: i, j, pnt

    ALLOCATE(pos(ndim))
    pos = pos0
    DO i = 1, ndim
       der(i) = 0.0_dp
       DO j = 1,-1,-2
          pos(i) = pos0(i) + j
          pnt    = point(pos,ndim,ngrid)
          der(i) = der(i) + REAL(j,KIND=dp)*(-fes(pnt))
       END DO
       pos(i) = pos0(i)
       der(i) = der(i) / (2.0_dp*dp_grid(i))
    END DO
    DEALLOCATE(pos)
  END FUNCTION derivative

END SUBROUTINE fes_min

! *****************************************************************************
!> \brief Finds path between two points (a) and (b)
!> \par History
!>      06.2009 created [tlaino]
!>      teodoro.laino .at. gmail.com 
!> \author Teodoro Laino
! *****************************************************************************
SUBROUTINE fes_path(fes, ndim, ngrid, dp_grid, x0, ndw, minima)
  USE kinds,                           ONLY: dp,&
                                             sp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndw
  INTEGER :: i,j,k,pnt,id,ntrials,nacc
  INTEGER, INTENT(in) :: ndim
  INTEGER, DIMENSION(:), POINTER :: pos, ngrid, pos0

  REAL(KIND=dp), DIMENSION(:), POINTER   :: xx, dx, rnd
  REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0
  REAL(KIND=dp)                          :: norm_dx, resto
  REAL(KIND=dp), DIMENSION(:,:), POINTER :: minima
  INTEGER, DIMENSION(ndim)               :: Dpos, ntrust
  LOGICAL                                :: do_save


END SUBROUTINE fes_path
