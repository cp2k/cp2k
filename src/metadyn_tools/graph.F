!-----------------------------------------------------------------------------!
!   FES: a fast and general program to map metadynamics on grids              !
!   Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013 !
!                 2014                                                        !
!                 Teodoro Laino                                               !
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
PROGRAM graph

  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE graph_methods,                   ONLY: fes_compute_low,&
                                             fes_cube_write,&
                                             fes_min,&
                                             fes_only_write,&
                                             fes_path,&
                                             fes_write
  USE graph_utils,                     ONLY: get_val_res,&
                                             mep_input_data_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_getarg,&
                                             m_iargc
  USE mathconstants,                   ONLY: pi

  IMPLICIT NONE

  CHARACTER(LEN=80)                        :: file, out1, out2, out3, wq_char, &
                                              path_file, out3_stride
  CHARACTER(LEN=default_string_length)     :: active_label, per_label
  INTEGER                                  :: argcount, coor, i, id, ip, &
                                              it, iw, ix, j, ncount, ndim, &
                                              ndw, nf, nfes, ngauss, nh, &
                                              nprd, nt, nt_p, nwr, p, q, &
                                              stat, unit_nr
  INTEGER, POINTER                         :: i_map(:), idw(:), ind(:), &
                                              inds(:), iperd(:), iprd(:), &
                                              ngrid(:), nn(:,:), nn_max(:), &
                                              tmp(:)
  LOGICAL                                  :: failure, fix, l_cp2k, l_cpmd, &
                                              l_dp, l_fes_int, l_fmin, &
                                              l_grid, l_math, l_orac, &
                                              l_pmin, lstride, l_popt, l_int,&
                                              l_cube
  REAL(KIND=dp)                            :: delta_s_glob, diff, dp2, dum, &
                                              eps_cut, sc, ss, x0w(3), &
                                              xfw(3)
  REAL(KIND=dp), POINTER :: delta_s(:,:), dp_cut(:), dp_grid(:), fes(:), &
    gauss(:,:), ss0(:,:), tmpr(:), ww(:), x0(:), xf(:)
  TYPE(mep_input_data_type)                :: mep_input_data


  failure=.FALSE.
  ! Initialize variables
  nprd    = 0
  ndim    = 1
  ndw     = 1
  nt_p    = 9999999
  eps_cut = 1e-6
  file    = 'HILLS'
  out1    = 'fes.dat'
  out2    = 'fes_int.dat'
  fix     = .FALSE.
  l_fes_int = .FALSE.
  lstride = .FALSE.
  l_grid  = .FALSE.
  l_dp    = .FALSE.
  l_orac  = .FALSE.
  l_cp2k  = .FALSE.
  l_cpmd  = .FALSE.
  l_math  = .FALSE.
  l_cube  = .FALSE.
  l_fmin  = .FALSE.
  l_pmin  = .FALSE.
  l_popt  = .FALSE.
  l_int   = .FALSE.
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
     WRITE(iw,*)'[-integrated_fes]    (When projecting the FES print the integrated value, '
     WRITE(iw,*)'                      rather then the minimum value (minimum value is default))'
     WRITE(iw,*)'[-orac]              (If energies are written in orac intern units)'
     WRITE(iw,*)'[-cp2k]              (Specify if a CP2K restart file is provided)'
     WRITE(iw,*)'[-cpmd]              (Specify if CPMD colvar_mtd and parvar_mtd are provided)'
     WRITE(iw,*)'                     (With CPMD you do not need to specify -file, parvar_mtd and'
     WRITE(iw,*)'                      colvar_mtd are expected to be present in the working directory)'
     WRITE(iw,*)'[-mathlab]           (File storing FES in Mathlab format. Default format Gnuplot)'
     WRITE(iw,*)'[-cube]              (File storing FES in GAUSSIAN CUBE format. Default format Gnuplot)'
     WRITE(iw,*)'[-find-minima]       (Tries to finds all minima in the computed FES)'
     WRITE(iw,*)'[-find-path]         (Finds MEP between all minima (found) in the computed FES)'
     WRITE(iw,*)'[-point-a]           (Specifies point (a) when using -find-path option)'
     WRITE(iw,*)'[-point-b]           (Specifies point (b) when using -find-path option)'
     WRITE(iw,*)'[-opt-path filename] (Optimize initial MEP of mep-nreplica points in the same format as mep.data)'
     WRITE(iw,*)'[-mep-kb]            (Specifies the value of the force constant for the MEP: default 0.1_dp)'
     WRITE(iw,*)'[-mep-nreplica]      (Specifies the number of replica points used in the MEP: default 8)'
     WRITE(iw,*)'[-mep-iter]          (Specifies the maximum number of iterations used in the MEP: default 10000)'
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
        READ(wq_char,*)out1
        ! we read only 1 filename. If none is specified we differentiate between fes.dat and fes_int.dat
        ! otherwise we use the one provided by the user
        out2=out1
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

     IF (INDEX(wq_char,'-integrated_fes').NE.0)THEN
        l_fes_int=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-orac').NE.0)THEN
        l_orac=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-cp2k').NE.0)THEN
        l_cp2k=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-cpmd').NE.0)THEN
        l_cpmd=.TRUE.
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

     IF (INDEX(wq_char,'-cube').NE.0)THEN
        l_cube=.TRUE.
     ENDIF

     IF (INDEX(wq_char,'-opt-path').NE.0)THEN
       l_popt=.TRUE.
       CALL M_GETARG(i+1,wq_char)
       READ(wq_char,*)path_file
     ENDIF

  END DO
  IF (COUNT((/l_orac,l_cp2k,l_cpmd/))/=1) &
       STOP "Error! You've to specify either ORAC, CP2K or CPMD!!"

  ! For CPMD move filename to colvar_mtd
  IF (l_cpmd) THEN
     file = "colvar_mtd"
  END IF

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
        l_dp  =.TRUE.
        l_grid=.FALSE.
        DO ix=1,ndim
           CALL M_GETARG(i+ix,wq_char)
           READ(wq_char,*)dp_grid(ix)
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
     ALLOCATE(mep_input_data%minima(ndw,2))
     mep_input_data%minima  =HUGE(0.0_dp)
     mep_input_data%max_iter=10000
     mep_input_data%kb=0.1_dp
     mep_input_data%nreplica=8
     ! Read for starting point (a) and (b)
     DO i=1,argcount
        CALL M_GETARG(i,wq_char)
        
        IF (INDEX(wq_char,'-point-a').NE.0)THEN
           DO id=1,ndw
              CALL M_GETARG(i+id,wq_char)
              READ(wq_char,*)mep_input_data%minima(id,1)
           ENDDO
        ENDIF

        IF (INDEX(wq_char,'-point-b').NE.0)THEN
           DO id=1,ndw
              CALL M_GETARG(i+id,wq_char)
              READ(wq_char,*)mep_input_data%minima(id,2)
           ENDDO
        ENDIF

        IF (INDEX(wq_char,'-mep-iter').NE.0)THEN
           CALL M_GETARG(i+1,wq_char)
           READ(wq_char,*)mep_input_data%max_iter
        ENDIF

        IF (INDEX(wq_char,'-mep-kb').NE.0)THEN
           CALL M_GETARG(i+1,wq_char)
           READ(wq_char,*)mep_input_data%kb
        ENDIF

        IF (INDEX(wq_char,'-mep-nreplica').NE.0)THEN
           CALL M_GETARG(i+1,wq_char)
           READ(wq_char,*)mep_input_data%nreplica
        ENDIF

     END DO
     IF (ANY(mep_input_data%minima==HUGE(0.0_dp))) &
          STOP "-find-path requires the specification of -point-a and -point-b !"
  ELSE
     ALLOCATE(mep_input_data%minima(0,0))
  END IF

! Read parameters for Path_optimization 
   IF (l_popt) THEN
        mep_input_data%nreplica=0
        mep_input_data%max_iter=10000
        mep_input_data%kb=0.1_dp

     DO i=1,argcount
        CALL M_GETARG(i,wq_char)

        IF (INDEX(wq_char,'-mep-kb').NE.0)THEN
           CALL M_GETARG(i+1,wq_char)
           READ(wq_char,*)mep_input_data%kb
        ENDIF

        IF (INDEX(wq_char,'-mep-iter').NE.0)THEN
           CALL M_GETARG(i+1,wq_char)
           READ(wq_char,*)mep_input_data%max_iter
        ENDIF

        IF (INDEX(wq_char,'-mep-nreplica').NE.0)THEN
           CALL M_GETARG(i+1,wq_char)
           READ(wq_char,*)mep_input_data%nreplica
        ENDIF
     END DO

        ALLOCATE(mep_input_data%minima(ndw,mep_input_data%nreplica))

        CALL open_file(unit_number=unit_nr,file_name=path_file,file_status='old')
        DO id=1, mep_input_data%nreplica
         READ(unit_nr,*) j, mep_input_data%minima(:,id)
        END DO
        CALL close_file(unit_nr)

        DO id=1, mep_input_data%nreplica
         WRITE(*,*)  mep_input_data%minima(:,id)
        END DO
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

  ! Revert the order so we can perform averages (when projecting FES) more
  ! efficiently
  i_map=ndim-i_map+1

  ! Tag the periodic COLVAR according the new internal order
  DO id=1,nprd
     iperd(i_map(iprd(id)))=1
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
  IF      (l_cp2k)            THEN
     CALL get_val_res(unit=10,section="&METADYN",keyword="NHILLS_START_VAL", i_val=nt)
     ! These sections may not necessarily be present.. if not the values will be HUGE and negative..
     CALL get_val_res(unit=10,section="&METADYN",keyword="HILL_TAIL_CUTOFF", r_val=sc)
     CALL get_val_res(unit=10,section="&METADYN",keyword="P_EXPONENT", i_val=p)
     CALL get_val_res(unit=10,section="&METADYN",keyword="Q_EXPONENT", i_val=q)
  ELSE IF (l_orac.OR.l_cpmd)  THEN
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

  IF      (l_cp2k) THEN
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_POS")
     DO i = 1, nt
        READ(10,*)(ss0(i_map(id),i),id=1,ndim)
     END DO
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_SCALE")
     DO i = 1, nt
        READ(10,*)(delta_s(i_map(id),i),id=1,ndim)
     END DO
     CALL get_val_res(unit=10,section="&METADYN",subsection="&SPAWNED_HILLS_HEIGHT")
     DO i = 1, nt
        READ(10,*)ww(i)
     END DO
  ELSE IF (l_orac) THEN
     DO i = 1, nt
        READ(10,*)dum,(ss0(i_map(id),i),id=1,ndim),(delta_s(i_map(id),i),id=1,ndim),ww(i)
     END DO
  ELSE IF (l_cpmd) THEN
     OPEN(11,file="parvar_mtd",status='old')
     DO i = 1, nt
        READ(10,*)dum,(ss0(i_map(id),i),id=1,ndim),(delta_s(id,i),id=1,ndim)
        READ(11,*)dum,dum,delta_s_glob,ww(i)
        delta_s(1:ndim,i)=delta_s_glob*delta_s(1:ndim,i)
     END DO
     CLOSE(11)
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
           IF (iperd(id)==1) THEN
              x0(id)=-pi
              xf(id)= pi
           ELSE
              x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
              xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
           END IF
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
     dp_grid=(xf-x0)/DBLE(ngrid-1)
  END IF
  
  WRITE(iw,'(70("*"))')
  WRITE(iw,'("FES|",T7,A,/)')"Parameters for FES:"
  WRITE(iw,'("FES|",T7,A15,5x,i7)')"NDIM         ::",ndim
  WRITE(iw,'("FES|",T7,A15,5x,i7)')"NWD          ::",ndw
  WRITE(iw,'("FES|",T7,A15,5x,i7)')"HILLS        ::",nt
  it=0
  DO i=ndim,1,-1
     it = it + 1
     per_label = ""
     active_label = "(NO MAPPED)"
     IF (iperd(i)/=0) per_label = "(PERIODIC)"
     IF (it<=ndw)     active_label = "(   MAPPED)"
     j = MINLOC((i_map-i)**2,1)
     WRITE(iw,'("FES|",T7,"COLVAR # ",i3," ::",5x,"(",f7.3," ,",f7.3,")",T48,A,T60,A)')&
          j,x0(i),xf(i),TRIM(per_label),TRIM(active_label)
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
              IF(sc > 0.0_dp) THEN
                 gauss(ip,i)= gauss(ip,i) * (1-(diff/sc*delta_s(i,it))**p)/(1-(diff/sc*delta_s(i,it))**q)
              END IF
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
        IF (l_fes_int) THEN
           out3=out2//"."
        ELSE
           out3=out1//"."
        END IF

        IF(ncount<10) THEN
           WRITE(out3_stride,'(A,i1)')TRIM(out3),ncount
        ELSEIF(ncount<100) THEN
           WRITE(out3_stride,'(A,i2)')TRIM(out3),ncount
        ELSE
           WRITE(out3_stride,'(A,i3)')TRIM(out3),ncount
        END IF
        OPEN(123,file=out3_stride)
        ind   = 1
        CALL fes_only_write(ndim, fes,  ind,  ndim, ngrid, dp_grid, ndw, l_fes_int)
        CLOSE(123)
     END IF

     nh=nh+nt_p
     nf=MIN(nh+nt_p-1,nt)
  END DO Stride
  DEALLOCATE(gauss)

  IF (l_fes_int) THEN
     out3=out2
  ELSE
     out3=out1
  END IF
 
  WRITE(iw,'("FES|",T7,A)') "Dumping FES structure in file: < "//TRIM(out3)//" >"
  OPEN(123,file=out3)
  IF (l_cube) THEN
     ind   = 1 
     CALL fes_cube_write(ndim, fes, ind, ndim, ngrid, dp_grid, x0, ndw, l_fes_int, file)
  ELSE
     ix=0
     IF (l_math) WRITE(123,'(10g12.5)')(ngrid(id),id=ndim,ndim-ndw+1,-1),ix
     ind   = 1   
     CALL fes_write(ndim, fes, ind, ndim, ngrid, dp_grid, x0, ndw, l_fes_int)
  END IF
  CLOSE(123)

  ! If requested find minima
  IF (l_fmin) CALL fes_min(fes, ndim, iperd, ngrid, dp_grid, x0, ndw)

  ! If requested find  or opt path
  IF ((l_pmin) .AND. (l_popt)) STOP "USE EITHER -find-path OR -opt-path"
  IF (l_pmin) l_int=.TRUE.
  IF (l_popt) l_int=.FALSE.

  IF ((l_pmin) .OR. (l_popt)) CALL fes_path(fes, ndim, ngrid, dp_grid, iperd, x0, ndw, mep_input_data, l_int)
  
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
  DEALLOCATE(mep_input_data%minima,stat=stat)
  IF (stat/=0) STOP "Deallocation Error"

  ! Terminate FES
  WRITE(iw,'(/,A,/)') "FES| NORMAL FES TERMINATION."

END PROGRAM graph
