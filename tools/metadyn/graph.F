!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
!-----------------------------------------------------------------------------!


!!****h* cp2k/fes [1.0] *
!!
!!   NAME
!!     fes
!!
!!   FUNCTION
!!     Program to Map on grid the hills spawned during a metadynamics run
!!
!!   NOTES
!!     Fast and improved version.
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [tlaino]
!!     teodoro.laino .at. gmail.com 
!!
!!
!!   SOURCE
!!     Please report any bug to the author
!!
!****************************************************************************

PROGRAM fes_compute

  USE kinds,                           ONLY: dp,&
                                             sp
  USE machine,                         ONLY: m_getlog,&
                                             m_hostnm,&
                                             m_getarg,&
                                             m_iargc
  
  IMPLICIT NONE
  
  INTERFACE
     SUBROUTINE get_val_res(unit, section, keyword, subsection, i_val)
       IMPLICIT NONE
       INTEGER, INTENT(IN)            :: unit
       CHARACTER (len=*)              :: section 
       CHARACTER (len=*),OPTIONAL     :: subsection 
       CHARACTER (len=*),OPTIONAL     :: keyword
       INTEGER, OPTIONAL, INTENT(OUT) :: i_val
     END SUBROUTINE get_val_res

     RECURSIVE SUBROUTINE fes_compute_low(idim, nn, fes, gauss, ind, ind0, nfes, ndim, ngauss, ngrid)
       USE kinds,                           ONLY: dp,&
                                                  sp       
       IMPLICIT NONE
       INTEGER, INTENT(in) :: idim, nfes, ngauss, ndim
       INTEGER, DIMENSION(:) :: nn, ind, ind0
       INTEGER, DIMENSION(:), POINTER :: ngrid
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes
       REAL(KIND=dp), DIMENSION(:,:), POINTER :: gauss
     END SUBROUTINE fes_compute_low

     RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw)
       USE kinds,                           ONLY: dp,&
                                                  sp
       IMPLICIT NONE
       INTEGER, INTENT(in) :: idim, ndim, ndw
       INTEGER, DIMENSION(:), POINTER :: pos
       INTEGER, DIMENSION(:), POINTER :: ngrid
       REAL(KIND=dp), DIMENSION(:), POINTER   :: fes, dp_grid, x0
     END SUBROUTINE fes_write

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
  END INTERFACE

  INTEGER :: i, ix ,id, it, ip, stat
  INTEGER :: ndim, ndw, nt_p, nt, nh, nf, nwr, ncount
  INTEGER :: argcount, ngauss, nfes, coor
  INTEGER, POINTER :: idw(:)
  INTEGER, POINTER :: i_map(:)
  INTEGER, POINTER :: ngrid(:), tmp(:)
  INTEGER, POINTER :: ind(:), inds(:), nn(:,:), nn_max(:)

  REAL(KIND=dp) :: dum, ss
  REAL(KIND=dp) :: x0w(3),xfw(3)
  REAL(KIND=dp), POINTER ::dp_cut(:)
  REAL(KIND=dp), POINTER :: ss0(:,:)
  REAL(KIND=dp), POINTER :: delta_s(:,:)
  REAL(KIND=dp), POINTER :: gauss(:,:),fes(:)
  REAL(KIND=dp), POINTER :: ww(:)
  REAL(KIND=dp), POINTER :: dp_grid(:),x0(:),xf(:), tmpr(:)
  REAL(KIND=dp) :: dp2,  eps_cut

  CHARACTER(LEN=80) :: wq_char, file, out

  LOGICAL :: fix, lstride
  LOGICAL :: l_grid, l_dp, l_orac, l_cp2k

  !
  ! Initialize variables
  !
  NDIM    = 2
  NDW     = 2
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


  argcount = m_iargc()
  IF(argcount==0)THEN
     WRITE(6,*)'USAGE:'
     WRITE(6,*)'graf  '
     WRITE(6,*)'[-ngrid 50   ]   (mesh dimension. DEFAULT :: 100)'
     WRITE(6,*)'[-ndim 3     ]   (number of collective variables NCV)'
     WRITE(6,*)'[-ndw 1 3    ]   (CVs for the free energy surface)'
     WRITE(6,*)'[-stride 10  ]   (how often the FES is written)'
     WRITE(6,*)'[-fix 1.1 ...]   (define the region for the FES)'
     WRITE(6,*)'                 (if omitted this is automatically calculated)'
     WRITE(6,*)'[-cutoff 2.]     (the hills are cutoffed at 2)'
     WRITE(6,*)'[-file filename]'
     WRITE(6,*)'[-out  filename]'
     WRITE(6,*)'[-orac ]         (if energies are written in orac intern units)'
     WRITE(6,*)'[-cp2k ]         (if CP2K restart file is provided)'
     WRITE(6,*)''
     WRITE(6,*)'DEFAULT OUTPUT: fes.dat'
     WRITE(6,*)''
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
        READ(wq_char,*)NDIM
     ENDIF

     IF (INDEX(wq_char,'-stride').NE.0)THEN
        CALL M_GETARG(i+1,wq_char)
        READ(wq_char,*)NT_P
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
  END DO
  IF (l_cp2k) l_orac = .FALSE.
  !
  ! Basic Allocation
  !
  ALLOCATE(ngrid(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE(dp_grid(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE(idw(ndw),stat=stat)
  IF (stat/=0) STOP "Allocation Error"

  DO i=1,ndw
     idw(i)=i
  END DO

  DO i=1,argcount
     CALL M_GETARG(i,wq_char)

     IF (INDEX(wq_char,'-ndw').NE.0)THEN
        NDW=0
        ndw_loop: DO ix=i+1,argcount
           CALL M_GETARG(ix,wq_char)
           IF(INDEX(wq_char,'-').EQ.0)THEN
              NDW=NDW+1
           ELSE
              EXIT ndw_loop
           ENDIF
        ENDDO ndw_loop
        DO id=1,NDW
           CALL M_GETARG(i+id,wq_char)
           READ(wq_char,*)idw(id)
        ENDDO
     ENDIF

     IF (INDEX(wq_char,'-ngrid').NE.0)THEN
        DO ix=1,ndim
           CALL M_GETARG(i+ix,wq_char)
           READ(wq_char,*)NGRID(ix)
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
        DO id=1,NDW
           CALL M_GETARG(i+2*(id-1)+1,wq_char)
           READ(wq_char,*)x0w(id)
           CALL M_GETARG(i+2*(id-1)+2,wq_char)
           READ(wq_char,*)xfw(id)
        ENDDO
     ENDIF
  ENDDO
  !
  !  defines the order of the collectiv var.: 
  !  first the "wanted" ones, then the others
  !  
  ALLOCATE(i_map(ndim),stat=stat)
  IF (stat/=0) STOP "Allocation Error" 

  DO id=1,NDW
     i_map(idw(id))=id
  ENDDO
  ix=NDW
  DO id=1,NDIM
     IF(i_map(id)==0)THEN
        ix=ix+1
        i_map(id)=ix
     ENDIF
  ENDDO
  i_map=ndim-i_map+1

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

  OPEN(10,file=file,status='old')
  IF (l_cp2k) THEN
     CALL get_val_res(unit=10,section="&METADYN",keyword="NHILLS_START_VAL", i_val=nt)
  ELSE
     OPEN(10,file=file,status='old')
     nt=0
     DO WHILE (.TRUE.)
        READ(10,*,END=100,ERR=100)dum
        nt=nt+1
     END DO
100  REWIND(10)
  END IF

  ALLOCATE( x0(NDIM)         , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( xf(NDIM)         , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( ss0(NDIM,NT)     , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( delta_s(NDIM,NT) , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( ww(NT)           , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( ind(NDIM)        , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( inds(NDIM)       , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( nn(NDIM,nt)      , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( nn_max(NDIM)     , stat=stat )
  IF (stat/=0) STOP "Allocation Error"
  ALLOCATE( dp_cut(NDIM)     , stat=stat )
  IF (stat/=0) STOP "Allocation Error"

  IF (l_cp2k) THEN
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_POS")
     DO i = 1, nt
        READ(10,*)(ss0(i_map(id),i),id=1,NDIM)
     END DO
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_AMPL")
     DO i = 1, nt
        READ(10,*)(delta_s(i_map(id),i),id=1,NDIM)
     END DO
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_HEIGHT")
     DO i = 1, nt
        READ(10,*)ww(i)
     END DO
  ELSE
     DO i = 1, nt
        READ(10,*)dum,(ss0(i_map(id),i),id=1,NDIM),&
                      (delta_s(i_map(id),i),id=1,NDIM),ww(i)
     END DO
  END IF
  CLOSE(10)

  IF(l_orac) ww = ww * 10000. / 4.187

  DO id=1,NDIM
     x0(id)=1000000
     xf(id)=-1000000
  ENDDO

  IF(fix) THEN
     DO it=1,NT
        DO id=1,NDIM-NDW
           x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
           xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
        ENDDO
     ENDDO
     it=0
     DO id=NDIM,NDIM-NDW+1,-1
        it=it+1
        x0(id)=x0w(it)
        xf(id)=xfw(it)
     ENDDO
  ELSE
     DO it=1,NT
        DO id=ndim,1,-1
           x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
           xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
        ENDDO
     ENDDO
  ENDIF
  
  IF(l_dp)THEN
     ALLOCATE(tmpr(ndim))
     tmpr=dp_grid
     DO i=1,ndim
        dp_grid(i_map(i))=tmpr(i)
     END DO
     DEALLOCATE(tmpr)
     ngrid=INT((xf-x0)/dp_grid)+1
  ELSE
     dp_grid=(xf-x0)/DBLE(NGRID-1)
  END IF
  
  WRITE(*,'(a15,5x,i7)'     )"NDIM         ::",NDIM
  WRITE(*,'(a15,5x,i7)'     )"NWD          ::",NDW
  WRITE(*,'(a15,5x,i7)'     )"HILLS        ::",NT
  it=ndim
  ip=1
  DO i=ndim-ndw+1,ndim
     WRITE(*,'(a9,i3,a3,5x,2f7.3)')"RANGE VAR",ip," ::",x0(it),xf(it)
     it=it-1
     ip=ip+1
  END DO
  WRITE(*,'(a15,5x,7i7)')     "NGRID        ::",(NGRID(id),id=ndim,ndim-ndw+1,-1)
  WRITE(*,'(a15,5x,5f7.3)')   "DX           ::",(dp_grid(id),id=ndim,ndim-ndw+1,-1)
  WRITE(*,'(a15,5x,10g10.5)' )"CUTOFF       ::",eps_cut

  nn_max = 0
  DO i = 1, nt
     dp_cut  = SQRT(LOG(ABS(ww(i))/eps_cut))*2.0_dp*delta_s(:,i) 
     nn(:,i) = INT(dp_cut/dp_grid)
     ww(i)   = ww(i)**(1.0_dp/DBLE(NDIM))
  END DO

  nn_max = MAXVAL(nn,DIM=2)
  ngauss = MAXVAL(nn_max) * 2 + 1
  nfes   = PRODUCT(ngrid)

  ALLOCATE(gauss(-MAXVAL(nn_max):MAXVAL(nn_max),NDIM))
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

  WRITE(*,'(a22)') "==>  Adding Hills  <=="
  stride : DO WHILE (nh <= nt)

     hills : DO it=nh,nf
        
        ind=INT((ss0(:,it)-x0)/dp_grid) + 1
        gauss=0.0_dp
        
        DO i=1,ndim
           coor = ind(i) - nn(i,it) - 1
           ss  = x0(i) + coor * dp_grid(i) - dp_grid(i)
           DO ip=-nn(i,it),nn(i,it)
              coor = coor + 1
              ss = ss + dp_grid(i)
              IF (coor .GT. ngrid(i)) CYCLE
              IF (coor .LT. 1) CYCLE
              dp2=((ss-ss0(i,it))/delta_s(i,it))**2
              gauss(ip,i)=ww(it)*EXP(-0.5_dp*dp2)
!              WRITE(*,*)i,ip, gauss(ip,i), ind, coor, ss, x0, ss0(i,it)
           END DO
        END DO
        inds = ind
        CALL fes_compute_low(NDIM,nn(:,it),fes,gauss,ind,inds,nfes,ndim,ngauss,ngrid)
        
        IF(.NOT. lstride .AND. MOD(it,nwr)==0)THEN
           WRITE(6,'(a13,i4,a2)') "Hills done ::",INT(10*ANINT(10.*it/nt))," %"
        ELSEIF(.NOT. lstride .AND. it==nt)THEN
           WRITE(6,'(a13,i4,a2)') "Hills done ::",INT(10*ANINT(10.*it/nt))," %"
        END IF

     END DO hills

     IF (lstride) THEN
        ncount = ncount+1
        WRITE(*,'(a13,i5," |-| Hills from ",i6," to",i6)') "Done frame ::",ncount,nh,nf
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
     
  END DO stride
  DEALLOCATE(gauss)

  WRITE(*,'(a22)') "==> Writing output <=="
  OPEN(123,file=out)
  ix=0
  WRITE(123,'(10g12.5)')(NGRID(id),id=ndim,ndim-ndw+1,-1),ix
  ind   = 1   
  CALL fes_write(ndim, fes,  ind, ndim, ngrid, dp_grid, x0, ndw)
  CLOSE(123)
  STOP "NORMAL END"
    
END PROGRAM fes_compute

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

RECURSIVE SUBROUTINE fes_compute_low(idim, nn, fes, gauss, ind, ind0, nfes, ndim, ngauss, ngrid)
  USE kinds,                           ONLY: dp,&
                                             sp       
  IMPLICIT NONE
  INTEGER :: i,pnt,j,k
  INTEGER, INTENT(in) :: idim, nfes, ngauss, ndim
  INTEGER, DIMENSION(:) :: nn, ind, ind0
  INTEGER, DIMENSION(:), POINTER :: pos , ll, ngrid
  
  REAL(KIND=dp), DIMENSION(:), POINTER   :: fes
  REAL(KIND=dp), DIMENSION(:,:), POINTER :: gauss
  REAL(KIND=dp) :: prod

  ALLOCATE(pos(ndim),ll(ndim))
  pos=ind
  k=nn(idim)

  DO i = -k, k
     pos(idim)=ind(idim)+i
     IF (pos(idim) .GT. ngrid(idim)) CYCLE
     IF (pos(idim) .LT. 1) CYCLE
     IF(idim/=1) THEN 
        CALL fes_compute_low(idim-1, nn, fes, gauss, pos, ind0, nfes, ndim, ngauss, ngrid)
     ELSE
        pnt=point(pos)
        prod=1.
        ll=pos-ind0
        DO j=1,ndim
           prod=prod*gauss(ll(j),j)
        END DO
        fes(pnt)=fes(pnt)+prod
     END IF
  END DO
  DEALLOCATE(pos)

  RETURN
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE fes_compute_low

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
        dimval = PRODUCT(NGRID(1:NDIM-NDW)) 
        WRITE(123,'(10f20.10)')(xx(id),id=ndim,ndim-ndw+1,-1),MINVAL(-fes(pnt:pnt+dimval-1))
     END IF     
  END DO
  DEALLOCATE(xx)
  RETURN
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE FES_WRITE

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
        dimval = PRODUCT(NGRID(1:NDIM-NDW)) !NGRID**(NDIM-NDW)
        WRITE(123,'(1f12.5)')MINVAL(-fes(pnt:pnt+dimval-1))
     END IF
     
  END DO
  RETURN
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE FES_ONLY_WRITE
