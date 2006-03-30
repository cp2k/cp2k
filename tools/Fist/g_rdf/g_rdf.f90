! Author : Teodoro Laino - SNS - ETHZ (Lugano)
! Date   : 26.11.2004. Released under LGPL license 2004.

PROGRAM main
  USE f90_unix
  USE f90_unix_dir
  USE f90_unix_env,                    ONLY: gethostname,&
                                             getlogin
  USE f90_unix_proc

  IMPLICIT NONE
  INTERFACE
     SUBROUTINE MakeSel(atom, residue, num, numres, element, ntap, sel, isel, verbose)
       IMPLICIT NONE
       CHARACTER (LEN=80)                       :: sel
       CHARACTER (LEN=5), POINTER, DIMENSION(:) :: atom 
       CHARACTER (LEN=3), POINTER, DIMENSION(:) :: residue
       CHARACTER (LEN=2), POINTER, DIMENSION(:) :: element
       INTEGER,           POINTER, DIMENSION(:) :: num
       INTEGER,           POINTER, DIMENSION(:) :: numres
       INTEGER,           POINTER, DIMENSION(:) :: isel
       INTEGER                                  :: ntap
       LOGICAL                                  :: verbose
     END SUBROUTINE MakeSel
  END INTERFACE
  
  INTEGER, PARAMETER  :: sp = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND ( 14, 200 )
  CHARACTER (LEN=80)  :: file_dcd, file_pdb, file_xsc, file_xyz, out_file
  CHARACTER (LEN=80)  :: sel1, sel2, diffusion
  LOGICAL             :: verbose
  REAL (KIND=sp), POINTER, DIMENSION(:,:,:)   :: coord
  CHARACTER (LEN=6), POINTER, DIMENSION(:)    :: name
  CHARACTER (LEN=5), POINTER, DIMENSION(:)    :: atom 
  CHARACTER (LEN=2), POINTER, DIMENSION(:)    :: element 
  CHARACTER (LEN=3), POINTER, DIMENSION(:)    :: residue 
  INTEGER, POINTER, DIMENSION(:)              :: num 
  INTEGER, POINTER, DIMENSION(:)              :: numres 
  INTEGER, POINTER, DIMENSION(:)              :: isel1, isel2
  REAL (KIND=dp), POINTER, DIMENSION(:)       :: gval, cell
  REAL(KIND=sp), POINTER, DIMENSION(:) ::  diff
  INTEGER        :: nframe, ntap, skip_frame, start_frame, end_frame, ref_frame, npoints
  REAL (KIND=dp) :: dr, dr_input, max_range

  NULLIFY( coord, name, atom, element, residue, num, numres, isel1, isel2, gval, cell, diff)

  CALL read_command_line(file_pdb, file_dcd, file_xsc, file_xyz, &
                         out_file, sel1, sel2, skip_frame,&
                         start_frame, end_frame, ref_frame, dr_input, diffusion,&
                         max_range, npoints, verbose)
  
  OPEN ( 9, file=file_pdb, status='old', form='formatted'  )
  IF (INDEX(file_dcd,"NULL") == 0) OPEN (10, file=file_dcd, status='old', form='unformatted')
  IF (INDEX(file_xyz,"NULL") == 0) OPEN (10, file=file_xyz, status='old', form='formatted')
  OPEN (11, file=file_xsc, status='old', form='formatted'  )
  
  CALL read_pdb( 9, name, atom, residue, num, numres, element, ntap)
  IF (INDEX(file_dcd,"NULL")  == 0) CALL read_dcd(coord, ntap, nframe,10)
  IF (INDEX(file_xyz,"NULL")  == 0) CALL read_xyz(coord, ntap, nframe,10)
  IF (INDEX(diffusion,"NULL") == 0) &
       CALL eval_diffusion(coord, atom, diffusion,start_frame, end_frame, skip_frame, ref_frame,&
       nframe, diff)
  CALL read_xsc(cell,11)
  ! Evaluate g(r)
  CALL makeSel(atom, residue, num, numres, element, ntap, sel1, isel1, verbose)
  CALL makeSel(atom, residue, num, numres, element, ntap, sel2, isel2, verbose)
  CALL evalG(coord, isel1, isel2, gval, cell, dr, nframe, &
       start_frame, end_frame, skip_frame, max_range, npoints)
  CALL DumpG(gval, file_pdb, out_file, sel1, sel2, dr)  

  DEALLOCATE(cell)
  DEALLOCATE(gval)
  DEALLOCATE(isel1, isel2)
  DEALLOCATE(coord)
  DEALLOCATE(name)
  DEALLOCATE(atom)
  DEALLOCATE(residue)
  DEALLOCATE(num)
  DEALLOCATE(numres)
  
  CLOSE( 9)
  CLOSE(10)
  CLOSE(11)
  
CONTAINS
  
  SUBROUTINE eval_diffusion(coord, atom, diffusion, start_frame, end_frame, skip_frame, ref_frame, nframe, diff)
    IMPLICIT NONE
    REAL(KIND=sp), POINTER, DIMENSION(:,:,:) :: coord
    REAL(KIND=sp), POINTER, DIMENSION(:) ::  diff
    CHARACTER (LEN=5), POINTER, DIMENSION(:) :: atom 
    CHARACTER (LEN=80)  :: diffusion
    INTEGER             :: skip_frame, start_frame, end_frame, ref_frame
    INTEGER             :: iframe, idim, j, nframe
    REAL(KIND=sp), POINTER, DIMENSION(:) :: x, y, z  

    IF (end_frame .LT. 0) end_frame = nframe
    idim = (end_frame - start_frame)/skip_frame + 1
    ALLOCATE(diff(idim))

    DO iframe = start_frame, end_frame, skip_frame
       x => coord(:,iframe,1)
       y => coord(:,iframe,2)
       z => coord(:,iframe,3)
       diff(iframe) = 0.0_sp
       DO j = 1, SIZE(atom)
          IF (TRIM(diffusion) == TRIM(atom(j))) &
               diff(iframe) = diff(iframe) + (coord(j,iframe,1)-coord(j,ref_frame,1))**2 +&
                                             (coord(j,iframe,2)-coord(j,ref_frame,2))**2 +&
                                             (coord(j,iframe,3)-coord(j,ref_frame,3))**2
       END DO
       diff(iframe) = SQRT(diff(iframe))
    END DO

    OPEN (12, file="diffusion", status="unknown", form="formatted")
    WRITE(12,'(A)')"# RMSD for atom ::"//TRIM(diffusion)
    DO j = 1, SIZE(diff)
       WRITE(12, '(f12.6)')diff(j)
    END DO
    CLOSE(12)
    DEALLOCATE(diff)
  END SUBROUTINE eval_diffusion

  SUBROUTINE print_help_banner
    IMPLICIT NONE
    
    WRITE(*,'(A)')"Usage: g_rdf -pdb filename -dcd filename -xsc filename"&
                            //" -sel selection [options]"
    WRITE(*,'(A)')"Options:"
    WRITE(*,'(2X,A,T20,A)')"-pdb","PDB file specifying types, resname"
    WRITE(*,'(2X,A,T20,A)')"-dcd","DCD file containing trajectory data (or XYZ)"
    WRITE(*,'(2X,A,T20,A)')"-xyz","XYZ file containing trajectory data (or DCD)"
    WRITE(*,'(2X,A,T20,A)')"-xsc","XSC file containing cell information"
    WRITE(*,'(2X,A,T20,A)')"-out","Optional. Specify the file name of the output data"
    WRITE(*,'(2X,A,T20,A)')"    ","Default is the pdb file's root name"
    WRITE(*,'(2X,A,T20,A)')"-sel","Specify selection for atom g(r)"
    WRITE(*,'(2X,A,T20,A)')"    ","2 selections needs to be specified"
    WRITE(*,'(2X,A,T20,A)')"-start_frame","Start Frame for evaluation of g(r)" 
    WRITE(*,'(2X,A,T20,A)')"-end_frame","End Frame for evaluation of g(r)" 
    WRITE(*,'(2X,A,T20,A)')"-skip_frame","Frames Stride" 
    WRITE(*,'(2X,A,T20,A)')"-diffusion","Evaluates the diffusion coefficient for the input TYPE" 
    WRITE(*,'(2X,A,T20,A)')"-max_range","Defines the max range on which the g(r) is evaluated"
    WRITE(*,'(2X,A,T20,A)')"-plotpoints ","Defines the number of points used in the range (see also -dr)"
    WRITE(*,'(2X,A,T20,A)')"-dr ","Defines the length of dr in the g(r) (see also -plotpoints)"
    WRITE(*,'(2X,A,T20,A)')"-verbose","Print verbose information during executions" 
    WRITE(*,'(2X,A,T20,A)')"    ","    "
    WRITE(*,'(2X,A,T20,A)')"Comments:",&
         "This program evaluates the g(r) given two atomic selections."
    STOP
  END SUBROUTINE print_help_banner

  SUBROUTINE read_command_line(file_pdb, file_dcd, file_xsc, file_xyz, out_file,&
                              sel1, sel2, skip_frame, start_frame, end_frame, ref_frame,&
                              dr_input, diffusion, max_range, npoints, verbose)
    IMPLICIT NONE
    CHARACTER(len=80)  :: file_dcd, file_pdb, file_xsc, out_file, file_xyz
    CHARACTER(len=80)  :: sel1, sel2, argval, diffusion
    INTEGER            :: skip_frame, start_frame, end_frame, ref_frame,narg, I, npoints
    LOGICAL            :: sel_set, verbose
    REAL(KIND=dp)      :: dr_input, max_range

    skip_frame  = 1
    start_frame = 1
    ref_frame   = -999
    max_range   = -999.0_dp
    npoints     =  999
    dr_input    =  0.1_dp
    end_frame   = -999
    verbose     = .FALSE.
    file_pdb    = "NULL"
    file_dcd    = "NULL"
    file_xsc    = "NULL"
    file_xyz    = "NULL"
    out_file    = "NULL"
    sel1        = "NULL"
    sel2        = "NULL"
    diffusion   = "NULL"
    sel_set     = .false.
    ! number of arguments...
    narg = iargc()
    IF (narg.LT.10) THEN
       CALL print_help_banner
    END IF
    DO i=1, narg
       CALL getarg(i, argval)
       SELECT CASE (TRIM(ARGVAL))
       CASE ("-dr")
          CALL getarg(i+1, argval)
          READ(argval,*)dr_input
       CASE ("-pdb")
          CALL getarg(i+1, file_pdb)
       CASE ("-dcd")
          CALL getarg(i+1, file_dcd)
       CASE ("-xyz")
          CALL getarg(i+1, file_xyz)
       CASE ("-xsc")
          CALL getarg(i+1, file_xsc)
       CASE ("-out")
          CALL getarg(i+1, out_file)
       CASE ("-sel")
          IF (sel_set) THEN
             CALL getarg(i+1,sel2)
          ELSE
             CALL getarg(i+1,sel1)
             sel_set = .TRUE.
          END IF
       CASE ("-start_frame")
          CALL getarg(i+1,argval)
          READ(argval,*)start_frame
       CASE ("-ref_frame")
          CALL getarg(i+1,argval)
          READ(argval,*)ref_frame
       CASE ("-end_frame")
          CALL getarg(i+1,argval)
          READ(argval,*)end_frame
       CASE ("-skip_frame")
          CALL getarg(i+1,argval)
          READ(argval,*)skip_frame
       CASE ("-diffusion")
          CALL getarg(i+1,diffusion)
       CASE ("-max_range")
          CALL getarg(i+1,argval)
          READ(argval,*)max_range
       CASE ("-plotpoints")
          CALL getarg(i+1,argval)
          READ(argval,*)npoints          
       CASE ("-verbose")
          verbose = .TRUE.
       CASE DEFAULT
          ! do nothing...
       END SELECT
    END DO

    IF (ref_frame < 0) ref_frame = start_frame
    IF (INDEX(file_pdb,"NULL") /= 0) CALL print_help_banner
    IF (INDEX(file_dcd,"NULL") /= 0 .AND. INDEX(file_xyz,"NULL") /= 0 ) CALL print_help_banner
    IF (INDEX(file_dcd,"NULL") == 0 .AND. INDEX(file_xyz,"NULL") == 0 ) CALL print_help_banner
    IF (INDEX(file_xsc,"NULL") /= 0) CALL print_help_banner
    IF (INDEX(    sel1,"NULL") /= 0) CALL print_help_banner
    IF (INDEX(    sel2,"NULL") /= 0) CALL print_help_banner
    
  END SUBROUTINE read_command_line
  

  SUBROUTINE DumpG(gval, file_pdb, out_file, sel1, sel2, dr)
    IMPLICIT NONE
    REAL(KIND=dp), POINTER, DIMENSION(:) :: gval
    REAL(KIND=dp)             :: dr
    CHARACTER (LEN=280)       :: gfile
    CHARACTER (LEN=80)        :: file_pdb, out_file, sel1, sel2, user, hostname
    CHARACTER (LEN=200)       :: cwd
    CHARACTER (LEN=80)        :: cal_date
    CHARACTER (LEN=10)        :: time
    CHARACTER (LEN=5)         :: zone
    CHARACTER (LEN=8)         :: date
    INTEGER, DIMENSION(8)     :: values
    INTEGER                   :: I, len, ipdbn

    CALL DATE_AND_TIME(date=date, time=time, zone=zone, values=values)
    cal_date=date//" "//time
    CALL getlogin(user,len)    
    CALL gethostname(hostname,len)
    CALL getcwd(cwd,len)
    IF (TRIM(out_file) == "NULL") THEN
       ipdbn = LEN_TRIM(file_pdb)
       DO WHILE ((file_pdb(ipdbn:ipdbn) /= "/").AND.(ipdbn > 1))
          ipdbn = ipdbn - 1
       END DO
       gfile = TRIM(cwd)//"/"//&
               TRIM(file_pdb(ipdbn:INDEX(file_pdb(ipdbn:),".")+ipdbn-1))//&
               "gr.dat"
    ELSE
       gfile = TRIM(out_file)
    END IF
    ! Dump info...
    IF (verbose) WRITE(*,*)" Dumping information on file :: ",gfile
    OPEN(12, file=gfile, status='unknown', form='formatted')
    WRITE(12,'(A)')"# Pair Correlation Function - g(r) - between :: "//&
         TRIM(sel1)//" and "//TRIM(sel2)//" ."
    WRITE(12,'(A)')"# Generated by :: "//TRIM(user)//" on :: "//TRIM(cal_date)//" ."
    WRITE(12,'(A)')"# Executed on  :: "//TRIM(hostname)//" ."
    WRITE(12,'(A)')"# Starting PDB file :: "//TRIM(file_pdb)//" ."
    WRITE(12,'(A)')"#        dr            g(r)"
    DO i = 1, SIZE(gval)
       WRITE(12,'(2F15.9)')REAL(i-1,dp)*dr,gval(i)
    END DO
    CLOSE(12)
  END SUBROUTINE DumpG

  SUBROUTINE EvalG(coord, isel1, isel2, gval, cell, dr, nframe,&
       start_frame, end_frame,  skip_frame, max_range, npoints)
    IMPLICIT NONE
    REAL(KIND=sp), POINTER, DIMENSION(:,:,:) :: coord    
    INTEGER, POINTER, DIMENSION(:) :: isel1, isel2
    REAL(KIND=dp), POINTER, DIMENSION(:) :: gval, cell
    REAL(KIND=sp), POINTER, DIMENSION(:) :: x, y, z  
    REAL(KIND=dp) :: gmax, dr, gnorm, pi, rho1, rho2, vol, r, max_range
    REAL(KIND=sp) :: iframe_count
    INTEGER :: np, nframe, iframe, skip_frame, i, start_frame, end_frame, npoints
    
    IF (end_frame .LT. 0) end_frame = nframe
    pi    = ATAN(1.0_dp)*4.0_dp
    gmax  = MINVAL(cell)/2.0_dp
    IF (max_range/= 999.0_dp) THEN
       gmax = max_range
    END IF
    IF (npoints /= 999) THEN
       dr = gmax/REAL(npoints,dp)
       np = npoints
    ELSE
       dr    = dr_input
       np    = FLOOR(gmax/dr)
       gmax  = dr * REAL(np,dp)
    END IF
    !
    vol   = PRODUCT(cell)
    rho1  = REAL(SIZE(isel1),dp)/vol
    rho2  = REAL(SIZE(isel2),dp)/vol
    gnorm = 4.0_dp*pi*dr*rho2
    
    ALLOCATE(gval(np))
    gval = 0.0_dp
    
    iframe_count = 0.0_dp
    DO iframe = start_frame, end_frame, skip_frame
       iframe_count = iframe_count + 1
       WRITE(*,'(A,I6)')" EVALG   :: Evaluating g(r) for frame number ::",iframe
       x => coord(:,iframe,1)
       y => coord(:,iframe,2)
       z => coord(:,iframe,3)
       CALL EvalG_low(x,y,z,cell,isel1,isel2, gval, dr, gmax)
    END DO
    ! At the end we need to normalize the g(r)...
    gval = gval / iframe_count
    DO i = 2, np
       r = REAL(i-1,dp)*dr
       gval(i) = gval(i) / (gnorm*r**2*rho1*vol)
    END DO
  END SUBROUTINE EvalG

  SUBROUTINE EvalG_low(x,y,z,cell,isel1,isel2,gval,dr,gmax)
    IMPLICIT NONE       
    INTEGER, POINTER, DIMENSION(:) :: isel1, isel2
    REAL(KIND=dp), POINTER, DIMENSION(:) :: gval, cell
    REAL(KIND=sp), POINTER, DIMENSION(:) :: x, y, z
    INTEGER :: i, j, ip, jp, ig
    REAL(KIND=dp) :: xd, xs, yd, ys, zd, zs, gmax, r, dr
    
    DO i = 1, SIZE(isel1)
       ip = isel1(i)
       DO j = 1, SIZE(isel2)
          jp = isel2(j)
          
          IF (ip == jp) CYCLE
          
          xd = x(ip) - x(jp)
          xs = xd / cell(1)
          xd = (xs - ANINT(xs,dp))*cell(1)
          
          yd = y(ip) - y(jp)
          ys = yd / cell(2)
          yd = (ys - ANINT(ys,dp))*cell(2)
          
          zd = z(ip) - z(jp)
          zs = zd / cell(3)
          zd = (zs - ANINT(zs,dp))*cell(3)
          
          r  = SQRT(xd**2+yd**2+zd**2)
          IF ( r .LT. gmax) THEN
             ig = FLOOR(r/dr)+1
             gval(ig) = gval(ig) + 1.0_dp
          END IF
       END DO
    END DO
  END SUBROUTINE EvalG_low
  
  SUBROUTINE read_xsc(cell, unit)
    IMPLICIT NONE
    REAL(KIND=dp), POINTER, DIMENSION(:) :: cell
    INTEGER             :: unit, idum
    CHARACTER (LEN=80)  :: line
    REAL      (KIND=dp) :: xdum
    
    ALLOCATE(cell(3))
    READ(unit,'(A)')line
    READ(unit,'(A)')line
    READ(unit,'(A)')line
    READ(line,*)idum, cell(1), xdum, xdum,&
                      xdum, cell(2), xdum,&
                      xdum, xdum, cell(3)
    WRITE(*,'(A,3F12.6)')" XSCINFO :: Cell Parameters :: ",cell

  END SUBROUTINE read_xsc

  SUBROUTINE read_dcd(coord, ntap, nframe, unit)
    IMPLICIT NONE
    INTEGER  :: ntap,unit,nframe
    INTEGER  :: i,m
    CHARACTER (LEN=4) :: car4
    INTEGER  :: nstart,nsanc,ntitle,namin,iframe
    INTEGER, DIMENSION(5) :: i5
    INTEGER, DIMENSION(9) :: i9
    CHARACTER (LEN=80):: car(10)
    REAL(KIND=dp)     :: delta
    REAL(KIND=sp), POINTER, DIMENSION(:,:,:) :: coord
    REAL(KIND=sp), POINTER, DIMENSION(:,:)   :: ax, ay, az
    
    REWIND(unit)
    READ (unit,END=100,ERR=200) car4, nframe, nstart, nsanc,i5,namin,delta,i9
    READ (unit,END=100,ERR=200) ntitle,(car(i),i=1,ntitle)
    WRITE(*,'(2X,A80)')(car(i),i=1,ntitle)
    READ (unit,END=100,ERR=200) ntap
    
    WRITE(*,'(A,I5)')" DCDINFO :: Number of Atoms  ::",ntap,&
                     " DCDINFO :: Number of Frames ::",nframe
    ALLOCATE( coord(ntap, nframe, 3) )
    ax => coord(:,:,1)
    ay => coord(:,:,2)
    az => coord(:,:,3)

    DO iframe = 1, nframe
       READ(unit,END=100,ERR=200)(ax(m,iframe),m=1,ntap)
       READ(unit,END=100,ERR=200)(ay(m,iframe),m=1,ntap)
       READ(unit,END=100,ERR=200)(az(m,iframe),m=1,ntap)
    END DO
    
100 RETURN
200 WRITE(*,*)"ERROR READING DCD FILE..."
    STOP
  END SUBROUTINE read_dcd

  SUBROUTINE read_xyz(coord, ntap, nframe, unit)
    IMPLICIT NONE
    INTEGER  :: ntap,unit,nframe
    INTEGER  :: i, idum, m
    INTEGER  :: iframe, ind
    CHARACTER (LEN=80):: car
    REAL(KIND=sp), POINTER, DIMENSION(:,:,:) :: coord
    REAL(KIND=sp), POINTER, DIMENSION(:,:)   :: ax, ay, az
    
    REWIND(unit)
    READ (unit,*,END=100,ERR=200) ntap
    REWIND(unit)
    iframe = 0
    DO WHILE (.TRUE.)
       READ (unit,*,END=99,ERR=200) idum
       IF (idum /= ntap) GOTO 200
       READ (unit,fmt='(A)',ERR=200)car
       DO i = 1, ntap
          READ (unit,fmt='(A)',ERR=200)car          
       END DO
       iframe = iframe + 1
       WRITE(*,'(A,I10)')" XYZINFO :: frame number ::",iframe 
    END DO
99  REWIND(unit)
    nframe = iframe
    WRITE(*,'(A,I5)')" XYZINFO :: Number of Atoms  ::",ntap,&
                     " XYZINFO :: Number of Frames ::",nframe
    ALLOCATE( coord(ntap, nframe, 3) )
    ax => coord(:,:,1)
    ay => coord(:,:,2)
    az => coord(:,:,3)

    DO iframe = 1, nframe
       READ (unit,*,END=100,ERR=200) ntap
       READ (unit,fmt='(A)',END=100,ERR=200) car
       DO i = 1, ntap
          READ(unit,fmt='(A)',END=100,ERR=200) car
          ind = INDEX(car," ")
          IF (ind == 1) THEN
             DO m = 2, LEN(car)
                IF (car(m:m) /= " ") EXIT
             END DO
             IF (m == LEN(car)+1) GOTO 200
             ind = m + INDEX(car(m+1:)," ")
          END IF
          IF (ind /= 1) READ(car(ind:),*)ax(i,iframe),ay(i,iframe),az(i,iframe)
       END DO
    END DO
    
100 RETURN
200 WRITE(*,*)"ERROR READING XYZ FILE..."
    STOP
  END SUBROUTINE read_xyz
  
  SUBROUTINE read_pdb(unit, name, atom, residue, num, numres, element, numatom)
    IMPLICIT NONE
    
    INTEGER :: jj,numatom,kk, unit, i
    CHARACTER (LEN=6)  :: str  
    CHARACTER (LEN=5)  :: myatom
    CHARACTER (LEN=3)  :: myresidue
    CHARACTER (LEN=80) :: line
    LOGICAL :: isatom  
    CHARACTER (LEN=6), POINTER, DIMENSION(:) :: name
    CHARACTER (LEN=5), POINTER, DIMENSION(:) :: atom 
    CHARACTER (LEN=2), POINTER, DIMENSION(:) :: element
    CHARACTER (LEN=3), POINTER, DIMENSION(:) :: residue
    INTEGER, POINTER, DIMENSION(:) :: num
    INTEGER, POINTER, DIMENSION(:) :: numres
    REAL(KIND=dp) :: xdum, ydum, zdum, occdum, chdum
    
11  FORMAT(a6,i5,1x,a5,a3,2x,i4,4x,3f8.3,2f6.2,4x,a2)
    
    numatom=0
    jj=0
    DO WHILE (.TRUE.)
       READ(unit,'(a80)',END=200)line
       READ(line,'(A6)')STR
       IF(INDEX(str,"ATOM").NE.0)THEN
          READ(line,'(A6,I5)')str,kk
          numatom=numatom+1
          jj=kk
       ELSEIF (INDEX(str,"TER").NE.0) THEN
          ! skip the line
       ELSEIF (INDEX(str,"END").NE.0) THEN
          EXIT
       ENDIF
    END DO
200 REWIND(unit)
    WRITE(*,'(A,I5)')" PDBINFO :: The number of atoms is ::",numatom
    
    ALLOCATE(name(numatom))
    ALLOCATE(atom(numatom))
    ALLOCATE(element(numatom))
    ALLOCATE(residue(numatom))
    ALLOCATE(num(numatom))
    ALLOCATE(numres(numatom))
    
    jj=0
    kk=0
    DO WHILE (.TRUE.)
       READ (unit,'(a80)',END=600) line
       isatom=(INDEX(line,"ATOM").NE.0)
       IF (INDEX(line,"TER").NE.0) THEN
          ! Do nothing
       ELSEIF (INDEX(line,"END").NE.0) THEN
          EXIT
       ELSEIF((kk.GT.0 .AND. kk.LE.numatom) .AND. .NOT.isatom)THEN
          WRITE(6,*)"wrong or inconsistent pdb file"
          WRITE(6,*)line
          STOP
       ENDIF
       IF(isatom)THEN
          kk=kk+1
          READ (line,11) name(kk),num(kk),atom(kk),  &
               residue(kk),numres(kk),xdum,ydum,zdum,&
               occdum,chdum
          ! modify white spaces..
          i = 1
          myatom = atom(kk)
          atom(kk) = "     "
          DO WHILE (myatom(i:i) == " ")
             i = i + 1
          END DO
          atom(kk) = myatom(i:5)
          i = 1
          myresidue = residue(kk)
          residue(kk) = "   "
          DO WHILE (myresidue(i:i) == " ")
             i = i + 1
          END DO
          residue(kk) = myresidue(i:3)
          i = 1
       ENDIF
    END DO
600 RETURN
  END SUBROUTINE read_pdb
  
END PROGRAM main

RECURSIVE SUBROUTINE makeSel(atom, residue, num, numres, element, ntap, sel, isel, verbose)
  IMPLICIT NONE

  TYPE array
     INTEGER, POINTER, DIMENSION(:) :: isel
     CHARACTER (LEN=81) :: myselarg
  END TYPE array
  
  TYPE p_array
     TYPE(array), POINTER         :: psel
  END TYPE p_array
  
  LOGICAL :: verbose
  CHARACTER (LEN=80) :: sel
  CHARACTER (LEN=1), DIMENSION(80) :: whitespaces
  CHARACTER (LEN=5), POINTER, DIMENSION(:) :: atom 
  CHARACTER (LEN=3), POINTER, DIMENSION(:) :: residue
  CHARACTER (LEN=2), POINTER, DIMENSION(:) :: element
  INTEGER, POINTER, DIMENSION(:) :: num
  INTEGER, POINTER, DIMENSION(:) :: numres
  LOGICAL ::  COR, CAND, CNOT, ANYLOGIC
  CHARACTER (LEN=81) :: mysel
  INTEGER, POINTER, DIMENSION(:) :: isel
  TYPE(p_array), DIMENSION(:), POINTER :: iselarg
  INTEGER :: icount, icount_OR, icount_AND, ii, ntap, iline
  INTEGER :: i, j, k, iselsize, isize, istart
  LOGICAL :: found
  INTEGER, SAVE :: itab
  DATA itab /0/
  
  NULLIFY(iselarg)
  COR   =.FALSE.
  CAND  =.FALSE.
  CNOT  =.FALSE.
  ANYLOGIC =.FALSE.
  
  mysel = " "//sel(1:80)
  itab  = itab + 5
  whitespaces = " "
  !
  ! Are there any logic prepositions ?
  !
  IF (verbose) WRITE(*,*)whitespaces(1:itab),"Entering with mysel ::",mysel
  IF (INDEX(mysel," OR ") /= 0 ) THEN
     COR = .TRUE.
     ! count them...
     icount_OR = 1
     iline     = 1
     DO WHILE (INDEX(mysel(iline:)," OR ") /= 0 )
        icount_OR = icount_OR + 1
        iline = iline + INDEX(mysel(iline:)," OR ") + 3
     END DO
     IF (verbose) WRITE(*,*)whitespaces(1:itab),"OR",icount_OR
     ALLOCATE(iselarg(icount_OR))
     iline    = 1
     iselsize = 0
     DO I = 1, icount_OR
        ALLOCATE(iselarg(i)%psel      )
        NULLIFY (iselarg(i)%psel%isel )
        IF (i.EQ.icount_OR) THEN
           iselarg(i)%psel%myselarg = TRIM(mysel(iline:))
        ELSE
           iselarg(i)%psel%myselarg = TRIM(mysel(iline:iline+INDEX(mysel(iline:)," OR ")-1))
           iline = iline + INDEX(mysel(iline:)," OR ") + 3
        END IF
        CALL MakeSel(atom, residue, num, numres, element, ntap,&
             iselarg(i)%psel%myselarg, iselarg(I)%psel%isel, verbose)
        iselsize = iselsize + SIZE(iselarg(I)%psel%isel)
        IF (verbose) WRITE(*,*)whitespaces(1:itab),&
             "Selection performed with array dimension :: ", SIZE(iselarg(I)%psel%isel)
        IF (verbose) WRITE(*,*)whitespaces(1:itab),"Partial Overall OR size dimension :: ",iselsize
     END DO
     
  ELSEIF (INDEX(mysel," AND ") /= 0 ) THEN
     CAND = .TRUE.
     ! count them...
     icount_AND = 1
     iline      = 1
     DO WHILE (INDEX(mysel(iline:)," AND ") /= 0 )
        icount_AND = icount_AND + 1
        iline = iline + INDEX(mysel(iline:)," AND ") + 4
     END DO
     IF(verbose) WRITE(*,*)whitespaces(1:itab),"AND",icount_AND
     ALLOCATE(iselarg(icount_AND))
     iline    = 1
     iselsize = 0
     DO I = 1, icount_AND
        ALLOCATE(iselarg(i)%psel)
        NULLIFY( iselarg(i)%psel%isel )
        IF (i.EQ.icount_AND) THEN
           iselarg(i)%psel%myselarg = TRIM(mysel(iline:))
        ELSE
           iselarg(i)%psel%myselarg = TRIM(mysel(iline:iline+INDEX(mysel(iline:)," AND ")-1))
           iline = iline + INDEX(mysel(iline:)," AND ") + 4
        END IF
        CALL MakeSel(atom, residue, num, numres, element, ntap,&
             iselarg(i)%psel%myselarg, iselarg(I)%psel%isel, verbose)
        IF (verbose) WRITE(*,*)whitespaces(1:itab),&
             "Selection performed with array dimension :: ", SIZE(iselarg(I)%psel%isel)
     END DO
     
  ELSEIF (INDEX(mysel," NOT ") /= 0) THEN
     CNOT = .TRUE.
     ALLOCATE(iselarg(1))
     iline    = 1
     iselsize = 0
     ALLOCATE(iselarg(1)%psel)
     NULLIFY( iselarg(1)%psel%isel )
     iline = iline + INDEX(mysel(iline:)," NOT ") + 4
     iselarg(1)%psel%myselarg = TRIM(mysel(iline:))
     
     CALL MakeSel(atom, residue, num, numres, element, ntap,&
          iselarg(1)%psel%myselarg, iselarg(1)%psel%isel, verbose)
     iselsize = ntap -  SIZE(iselarg(1)%psel%isel)
     IF (verbose) WRITE(*,*)whitespaces(1:itab),&
          "Selection performed with array dimension :: ", SIZE(iselarg(1)%psel%isel)
     IF (verbose) WRITE(*,*)whitespaces(1:itab),"Definitive Overall NOT size dimension :: ",iselsize
  END IF
  
  IF (COR.OR.CAND.OR.CNOT) ANYLOGIC=.TRUE.
  
  IF (ANYLOGIC) THEN
     IF (COR) THEN
        !
        ! With OR operator we unite the contributions
        !
        ALLOCATE(isel(iselsize))
        isize = 1
        DO i = 1, icount_OR
           isel( isize:isize+SIZE(iselarg(i)%psel%isel)-1)=iselarg(i)%psel%isel
           isize = isize + SIZE(iselarg(i)%psel%isel)
           DEALLOCATE(iselarg(i)%psel%isel)
        END DO
        DEALLOCATE(iselarg)
        IF (verbose) WRITE(*,*)whitespaces(1:itab),"gone out from OR.."
     ELSE IF (CAND) THEN
        !
        ! With AND operator we intersect the contributions
        !        
        DO k = 2, icount_AND
           istart = 1
           CALL Sort(iselarg(k-1)%psel%isel, SIZE(iselarg(k-1)%psel%isel))
           CALL Sort(iselarg(k  )%psel%isel, SIZE(iselarg(k  )%psel%isel))
           DO i = 1, SIZE(iselarg(k)%psel%isel)
              found = .FALSE.
              DO j = istart, SIZE(iselarg(k-1)%psel%isel)
                 IF ( iselarg(k-1)%psel%isel(j) .GT. &
                      iselarg(k  )%psel%isel(i)      ) THEN
                    istart = j
                    EXIT
                 END IF
                 IF ( iselarg(k  )%psel%isel(i) .EQ. &
                      iselarg(k-1)%psel%isel(j)      ) THEN                    
                    found = .TRUE.
                    EXIT
                 END IF
              END DO
              IF (.NOT.found) iselarg(k  )%psel%isel(i) = 999999999
           END DO
           DEALLOCATE(iselarg(k-1)%psel%isel)
        END DO
        iselsize = 0
        DO i = 1, SIZE(iselarg(icount_AND)%psel%isel)
           IF (iselarg(icount_AND)%psel%isel(i).NE.999999999) iselsize = iselsize + 1
        END DO
        ALLOCATE(isel(iselsize))
        iselsize = 0
        DO i = 1, SIZE(iselarg(icount_AND)%psel%isel)
           IF (iselarg(icount_AND)%psel%isel(i).NE.999999999) THEN
              iselsize = iselsize + 1
              isel(iselsize) = iselarg(icount_AND)%psel%isel(i)
           END IF
        END DO
        DEALLOCATE(iselarg(icount_AND)%psel%isel)
        DEALLOCATE(iselarg)
        IF (verbose) WRITE(*,*)whitespaces(1:itab),"uscito da AND.."
     ELSE
        !
        ! With NOT we get the complement of the contribution
        !
        ALLOCATE(isel(iselsize))
        icount = 0
        External_loop: DO ii = 1, ntap
           DO j = 1, SIZE(iselarg(1)%psel%isel)
              IF (ii.EQ.iselarg(1)%psel%isel(j)) CYCLE External_loop
           END DO
           icount = icount + 1
           isel(icount) = ii
        END DO External_loop
        DEALLOCATE(iselarg(1)%psel%isel)
        DEALLOCATE(iselarg)        
        IF (verbose) WRITE(*,*)whitespaces(1:itab),"uscito da NOT.."
     END IF
     IF (verbose) WRITE(*,*)whitespaces(1:itab),"Array selection dimension :: ",SIZE(isel)
  ELSE
     IF (verbose) WRITE(*,*)whitespaces(1:itab),"Making the elementary selection.."
     CALL MakeSel_Low(mysel, isel)
  END IF
  itab = itab -5
  IF (itab == 0) & ! Exiting completely from MakeSel... Punch out some results...
       WRITE(*,'(A,I6)')" MAKESEL :: Selected atom :: "//TRIM(sel)//&
       ". Total number of selected atoms ::",SIZE(isel)
  RETURN

CONTAINS

  RECURSIVE SUBROUTINE MakeSel_Low(mysel, isel)
    CHARACTER (LEN=*)  :: mysel
    CHARACTER (LEN=80) :: myname, LBLID, LBIND
    INTEGER, POINTER, DIMENSION(:) :: isel, myisel, myisela, myiselb, myiselc
    INTEGER :: ii_start, ii_end, res_id, ii, icount, i_from, i_to, i_end
    INTEGER :: start_loop, end_loop, ind_star
    LOGICAL :: CNAME, CALLL, CRES, CNRES, CID

    CALLL =.FALSE.
    CNAME =.FALSE.
    CRES  =.FALSE.
    CID   =.FALSE.
    CNRES =.FALSE.
    NULLIFY(myisel, myisela, myiselb, myiselc)

    IF (INDEX(mysel," NAME") /= 0 ) THEN
       CNAME = .TRUE.
       LBLID = " NAME"
    END IF
    
    IF (INDEX(mysel," ALL") /= 0 ) THEN
       CALLL = .TRUE.
       LBLID = " ALL"
    END IF
    
    IF (INDEX(mysel," RESNAME") /= 0) THEN
       CRES  = .TRUE.
       LBLID = " RESNAME"
    END IF
    
    IF (INDEX(mysel," RESIDUE") /= 0) THEN
       CNRES = .TRUE.
       LBLID = " RESIDUE"
    END IF
    
    IF (INDEX(mysel," INDEX") /= 0) THEN
       CID = .TRUE.
       LBLID = " INDEX"
    END IF
    
    icount = 0
    ! Identify the selection type
    IF ( CNAME.OR.CRES.OR.CNRES.OR.CID ) THEN
       ii = INDEX(mysel,TRIM(LBLID)) + LEN_TRIM(LBLID)
       IF (verbose) WRITE(*,*)mysel
       DO WHILE (mysel(ii:ii) == " ")
          ii = ii+1
       END DO
       ii_start = ii
       ii_end = LEN_TRIM(mysel(ii:))+ii-1
       myname = mysel(ii_start:ii_end)

       IF (INDEX(myname,"FROM") /= 0) THEN
          IF ( CNAME.OR.CRES.OR.INDEX(myname," TO") == 0) THEN
             WRITE(*,*)"MAKESEL :: FATAL ERROR IN SELECTION ::"//TRIM(mysel)
             WRITE(*,*)"MAKESEL :: (1) THE STATEMENT from .. to .. IS ONLY ALLOWED WITH :: INDEX AND RESIDUE!"
             WRITE(*,*)"MAKESEL :: (2) from STATEMENT DOESN'T MATCH to...  "
             WRITE(*,*)"MAKESEL :: INTERRUPTING SELECTION EVALUATION..."
             STOP
          END IF
          i_from = INDEX(myname,"FROM")
          i_to   = INDEX(myname(i_from:)," TO") + i_from - 1 
          i_end  = i_to + LEN_TRIM(" TO")
          DO WHILE (mysel(i_end:i_end) == " ")
             i_end = i_end + 1
          END DO
          DO WHILE (mysel(i_end:i_end) /= " ")
             i_end = i_end + 1
          END DO
          
          READ(myname(i_from+LEN_TRIM("FROM"):i_to),*)start_loop
          READ(myname(i_to  +LEN_TRIM(" TO") :i_end),*)end_loop
          DO ii = start_loop, end_loop
             WRITE(LBIND,*)ii
             CALL MakeSel_Low(TRIM(LBLID)//" "//TRIM(LBIND), myisela)
             IF (ASSOCIATED(myiselb)) THEN
                ALLOCATE(myiselc(SIZE(myiselb)))
                myiselc = myiselb
                DEALLOCATE(myiselb)
                ALLOCATE(myiselb(SIZE(myisela)+SIZE(myiselc)))
                myiselb(1:size(myisela)) = myisela
                myiselb(SIZE(myisela)+1:SIZE(myisela)+size(myiselc)) = myiselc
                DEALLOCATE(myiselc)
                DEALLOCATE(myisela)
             ELSE
                ALLOCATE(myiselb(SIZE(myisela)))
                myiselb(1:SIZE(myisela)) = myisela
                DEALLOCATE(myisela)
             END IF
          END DO
          
          isel => myiselb
          IF (.NOT.empty(myname(1:i_from-1)//" "//myname(i_end:LEN_TRIM(myname(i_end:))))) &
               CALL MakeSel_Low(TRIM(LBLID)//" "//myname(1:i_from-1)//&
               " "//myname(i_end:LEN_TRIM(myname(i_end:))), myisel)
       ELSE
          ii = LEN_TRIM(myname)
          DO WHILE (myname(ii:ii) /= " " .AND. ii > 1)
             ii = ii-1
          END DO

          IF (ii .GT. 1) THEN
             IF (verbose) WRITE(*,*)whitespaces(1:itab),&
                  "MAKESEL :: RECURSIVE CALL TO MAKESEL WITH SELECTION ::"&
                  //TRIM(TRIM(LBLID)//" "//myname(1:ii))//" ."          
             CALL MakeSel_Low(TRIM(LBLID)//" "//myname(1:ii), myisel)
             myname = myname(ii+1:LEN_TRIM(myname))
             IF (verbose) WRITE(*,*)whitespaces(1:itab),&
                  "MAKESEL :: EXITING FROM RECURSIVE CALL... DIMENSION OF SELECTION ::",SIZE(myisel)," ."
             IF (verbose) WRITE(*,*)whitespaces(1:itab),&
                  "MAKESEL :: CONTINUING WITH SELECTION :: "//TRIM(myname)// " ."
          END IF
          ! Count the number of atoms
          icount = 0
          ind_star = INDEX(TRIM(myname),"*")
          IF     (CNAME) THEN
             IF (ind_star /= 0) THEN
                DO ii = 1, ntap
                   IF (INDEX(TRIM(myname(1:ind_star)),TRIM(atom(ii))) == 1)  icount=icount+1
                END DO                
             ELSE
                DO ii = 1, ntap
                   IF (TRIM(myname) == TRIM(atom(ii)))    icount=icount+1
                END DO
             END IF
          ELSEIF (CRES ) THEN
             IF (ind_star /= 0) THEN
                IF (INDEX(TRIM(myname(1:ind_star)),TRIM(residue(ii))) == 1) icount=icount+1
             ELSE             
                DO ii = 1, ntap
                   IF (TRIM(myname) == TRIM(residue(ii))) icount=icount+1
                END DO
             END IF
          ELSEIF (CNRES) THEN
             icount = 1
             READ(myname,*)res_id
          ELSEIF (CID  ) THEN
             icount = 1          
          END IF
          
          IF (ASSOCIATED(myisel)) THEN
             IF (verbose) WRITE(*,*)whitespaces(1:itab),&
                  "icount  ::",icount," size(myisel) ::",SIZE(myisel)
             icount = icount + SIZE(myisel)
          ELSE
             IF (verbose) WRITE(*,*)whitespaces(1:itab),&
                  "icount  ::",icount
          END IF

          ALLOCATE(isel(icount))
          ! Make real selection...
          icount = 0
          IF     (CNAME) THEN
             IF (ind_star /= 0) THEN
                DO ii = 1, ntap
                   IF (INDEX(TRIM(myname(1:ind_star)),TRIM(atom(ii))) == 1) THEN
                      icount = icount + 1
                      isel(icount) = ii
                   END IF
                END DO
             ELSE             
                DO ii = 1, ntap
                   IF (TRIM(myname) == TRIM(atom(ii))) THEN
                      icount = icount + 1
                      isel(icount) = ii
                   END IF
                END DO
             END IF
          ELSEIF (CRES ) THEN
             IF (ind_star /= 0) THEN
                DO ii = 1, ntap
                   IF (INDEX(TRIM(myname(1:ind_star)),TRIM(residue(ii))) == 1) THEN
                      icount = icount + 1
                      isel(icount) = ii
                   END IF
                END DO
             ELSE
                DO ii = 1, ntap
                   IF (TRIM(myname) == TRIM(residue(ii))) THEN
                      icount = icount + 1
                      isel(icount) = ii
                   END IF
                END DO
             END IF
          ELSEIF (CNRES) THEN
             CALL  Cres_Sel(icount, res_id, residue, myname, ii_start, ii_end, isel)
          ELSEIF (CID  ) THEN
             icount = 1
             READ(myname,*)isel(1)
          END IF
       END IF

       IF (ASSOCIATED(myisel)) THEN
          isel(icount+1:) = myisel
          DEALLOCATE(myisel)
       END IF
       
       IF     (CNAME.OR.CID) THEN
          IF (verbose) WRITE(*,*)whitespaces(1:itab),&
               "MAKESEL :: Selected atom :: "//TRIM(myname)//". Total number of selected atoms ::",SIZE(isel)
       ELSEIF (CRES.OR.CNRES) THEN
          IF (verbose) WRITE(*,*)whitespaces(1:itab),&
               "MAKESEL :: Selected residue :: "//TRIM(myname)//". Total number of selected atoms ::",SIZE(isel)
       END IF           
    END IF
        
    IF (CALLL) THEN
       ALLOCATE(isel(ntap))
       DO ii = 1, ntap
          isel(ii) = ii
       END DO
       IF (verbose) WRITE(*,*)whitespaces(1:itab),&
            "MAKESEL :: Selecting all atoms    . Total number of selected atoms ::",SIZE(isel)
    END IF
    
  END SUBROUTINE MakeSel_Low
   

  SUBROUTINE Cres_Sel(icount, res_id, residue, myname, ii_start, ii_end, isel)
    IMPLICIT NONE
    INTEGER :: icount, res_id, ii_start, ii_end
    CHARACTER (LEN=80) :: myname
    INTEGER :: ii
    INTEGER, POINTER, DIMENSION(:) :: isel
    CHARACTER (LEN=3), POINTER, DIMENSION(:) :: residue

    icount = 1
    ii     = 1
    myname = TRIM(residue(ii))
    DO WHILE (icount /= res_id)
       ii = ii + 1
       IF (myname /= residue (ii)) icount = icount + 1
    END DO
    icount = 0
    ii_start = ii
    myname = TRIM(residue(ii_start))       
    DO WHILE (TRIM(myname) == residue(ii))
       icount = icount + 1
       ii = ii + 1
    END DO
    ii_end = ii - 1
    ALLOCATE(isel(ii_end-ii_start+1))
    DO ii = ii_start, ii_end
       isel(ii-ii_start+1) = ii
    END DO
   
  END SUBROUTINE Cres_Sel

  LOGICAL FUNCTION EMPTY(string) RESULT(IS_EMPTY)
    IMPLICIT NONE
    INTEGER :: i
    CHARACTER (LEN=*) :: string

    IS_EMPTY = .TRUE.
    DO i = 1, LEN_TRIM(string)
       IF (string(i:i) /= " ") THEN
          is_empty = .FALSE.
          EXIT
       END IF
    END DO
  END FUNCTION EMPTY

  INTEGER FUNCTION FindMinimum(x, Start, END)
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN)               :: Start, END
    INTEGER                           :: Minimum, Location, I
    
    Minimum = x(Start)
    Location = Start
    DO I = Start + 1, END
       IF (x(I) < Minimum) THEN
          Minimum = x(I)
          Location = I
       END IF
    END DO
    FindMinimum = Location
  END FUNCTION FindMinimum
  
  SUBROUTINE Swap(A, B)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: A, B
    INTEGER                :: Tmp
    
    Tmp = A
    A = B
    B = Tmp
  END SUBROUTINE Swap
  
  SUBROUTINE Sort(x, Size)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                  :: Size
    INTEGER, DIMENSION(:), INTENT(INOUT) :: x
    INTEGER                              :: I, Location
    
    DO I = 1, Size - 1
       Location = FindMinimum(x, i, Size)
       CALL Swap(x(i), x(Location))
    END DO
  END SUBROUTINE Sort
  
END SUBROUTINE makeSel
