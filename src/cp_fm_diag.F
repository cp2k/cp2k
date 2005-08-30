!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_fm_diag [1.0] *
!!
!!   NAME
!!     cp_fm_diag
!!
!!   FUNCTION
!!     used for collecting some of the diagonalization shemes available for cp_fm_type
!!     cp_fm_power also moved here as it is very related
!!
!!   NOTES
!!     first version : most routines imported 
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
! ************************************************************************************

MODULE cp_fm_diag
  USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                             cp_blacs_gridinit
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_syrk
  USE cp_fm_types,                     ONLY: cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_unit_nr,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level, &
                                             cp_logger_get_default_unit_nr
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_diag'

!****************************************************************************

PUBLIC :: cp_fm_syevd, &
          cp_fm_syevx, &
          cp_fm_power, &
          cp_fm_block_jacobi_classic  ! more jacobi variants are present in the file but currently not used

CONTAINS

!*************************************************************************
!computes all eigenvalues and vectors of a real symmetric matrix
!should be quite a bit faster than syevx for that case
!especially in parallel with thightly clustered evals
!needs more workspace in the worst case, but much better distributed
!*************************************************************************
  SUBROUTINE cp_fm_syevd(matrix,eigenvectors,eigenvalues,error)

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syevd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, istat, liwork, &
                                                lwork, mypcol, myprow, n
    INTEGER, DIMENSION(9)                    :: descm, descv
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
    TYPE(cp_blacs_env_type), POINTER         :: context

!!FM    integer :: muno1,muno2

    CALL timeset("cp_fm_syevd","I","",handle)
    failure=.FALSE.

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = matrix%matrix_struct%nrow_global

    m => matrix%local_data

#if defined(__SCALAPACK)
    descm(:) = matrix%matrix_struct%descriptor(:)
    v => eigenvectors%local_data
    descv(:) = eigenvectors%matrix_struct%descriptor(:)

    liwork=7*n+8*context%num_pe(2)+2
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("cp_fm_syevd","iwork")
! work space query
    lwork = -1
    ALLOCATE(work(1),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("cp_fm_syevd","work")
    CALL PDSYEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,iwork(1),liwork,info)
    ! look here for a PDORMTR warning :-)
    ! this routine seems to need more workspace than reported
    ! (? lapack bug ...)
    ! arbitrary additional memory  ... we give 100000 more words (it seems to depend on the block size used)
    lwork = NINT(work(1)+100000)
    DEALLOCATE(work)
    ALLOCATE(work(lwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("cp_fm_syevd","lwork")
    CALL PDSYEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,iwork(1),liwork,info)
#else
    ! for n=4 with my lapack (linux) I need to add 50 to the value of lwork. 5000 seems to work in IRIX
    ! unfortunately query is not supported. I did not test n dependence[fawzi]
    lwork=1+6*n+2*n**2+5000
    liwork=5*n+3
    ALLOCATE(work(lwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("cp_fm_syevd","work")
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("cp_fm_syevd","iwork")
    CALL DSYEVD('V','U', n, m(1,1), n, eigenvalues(1), work(1), lwork, iwork(1), liwork, info)
    CALL cp_fm_to_fm(matrix,eigenvectors)
#endif

    IF (info.NE.0) CALL stop_program("cp_fm_syevd","unable to diagonalize matrix")

    DEALLOCATE(iwork)
    DEALLOCATE(work)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_syevd

! *****************************************************************************

  SUBROUTINE cp_fm_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,error)

!   Purpose: Diagonalise the symmetric n by n matrix using the LAPACK library.
!            if bmatrix present, use general solver, will overwrite bmatrix
!            with its cholesky decom.

!   History: - Creation (06.06.2001, Matthias Krack)
!   History: - Mod      (05.2001, Joost VandeVondele)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER          :: eigenvectors,matrix
    REAL(KIND = dp), INTENT(IN)                      :: work_syevx
    INTEGER, INTENT(IN)                       :: neig
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT)       :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_fm_syevx",&
      routineP=moduleN//':'//routineN
    REAL(KIND = dp), PARAMETER :: orfac = -1.0_dp,&
                           vl = 0.0_dp,&
                           vu = 0.0_dp

!   *** Local variables ***

    REAL(KIND = dp) :: abstol
    INTEGER  :: handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,&
                nn,np0,npcol,npe,nprow,nq0,nz,output_unit,itype
    LOGICAL  :: ionode

    INTEGER, DIMENSION(9) :: desca,descz,descb

    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: gap,w,work
    INTEGER, DIMENSION(:), ALLOCATABLE  :: iclustr,ifail,iwork
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: a,z,b
    TYPE(cp_blacs_env_type), POINTER    :: context
    TYPE(cp_logger_type), POINTER       :: logger

    REAL(KIND = dp), EXTERNAL :: dlamch

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: iceil,numroc
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif

!   ---------------------------------------------------------------------------

    IF (neig==0) RETURN

    CALL timeset("cp_fm_syevx","I","",handle)

    abstol = 2.0_dp*dlamch("S")

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source

    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    n = matrix%matrix_struct%nrow_global

    ALLOCATE (w(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w",n*dp_size)
    w(:) = 0.0_dp

    eigenvalues(:) = 0.0_dp
#if defined(__SCALAPACK)

    IF (matrix%matrix_struct%nrow_block /= matrix%matrix_struct%ncol_block) THEN
       CALL stop_program(routineP,"Invalid blocksize (no square blocks)")
    END IF

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)
    z => eigenvectors%local_data
    descz(:) = eigenvectors%matrix_struct%descriptor(:)

!   *** Get the optimal work storage size ***

    npe = nprow*npcol
    nb = matrix%matrix_struct%nrow_block
    nn = MAX(n,nb,2)
    np0 = numroc(nn,nb,0,0,nprow)
    nq0 = MAX(numroc(nn,nb,0,0,npcol),nb)

    lwork = 5*n + MAX(5*nn,np0*nq0) + iceil(neig,npe)*nn + 2*nb*nb +&
            INT(work_syevx*REAL((neig - 1)*n,dp)) !!!! allocates a full matrix on every CPU !!!!!
    liwork = 6*MAX(N,npe+1,4)

    ALLOCATE (gap(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"gap",npe*dp_size)
    gap = 0.0_dp
    ALLOCATE (iclustr(2*npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iclustr",2*npe*int_size)
    iclustr = 0
    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail",n*int_size)
    ifail = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork",liwork*int_size)
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",lwork*dp_size)

    CALL pdsyevx("V","I","U",n,a(1,1),1,1,desca,vl,vu,1,neig,abstol,m,nz,w(1),orfac,&
                 z(1,1),1,1,descz,work(1),lwork,iwork(1),liwork,ifail(1),iclustr(1),gap,info)


!   *** Error handling ***

    IF (info /= 0) THEN
      IF (ionode) THEN
         output_unit = cp_logger_get_unit_nr(logger,&
              cp_failure_level,routineP,local=.FALSE.)
        WRITE (unit=output_unit,FMT="(/,(T3,A,T12,1X,I10))")&
          "info    = ",info,&
          "lwork   = ",lwork,&
          "liwork  = ",liwork,&
          "nz      = ",nz
        IF (info > 0) THEN
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
            "ifail   = ",ifail
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
            "iclustr = ",iclustr
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,E10.3)))")&
            "gap     = ",gap
        END IF
      END IF
      CALL stop_program(routineP,"Error in pdsyevx")
    END IF

!   *** Release work storage ***

    DEALLOCATE (gap,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"gap")
    DEALLOCATE (iclustr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iclustr")
    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work")

#else

    a => matrix%local_data
    z => eigenvectors%local_data

!   *** Get the optimal work storage size ***

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail",n*int_size)
    ifail = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork",liwork*int_size)
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",lwork*dp_size)
    info = 0

    CALL dsyevx("V","I","U",n,a(1,1),n,vl,vu,1,neig,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)

!   *** Error handling ***

    IF (info /= 0) THEN
      output_unit = cp_logger_get_unit_nr(logger,&
                    cp_failure_level,routineP,local=.FALSE.)
      WRITE (unit=output_unit,FMT="(/,(T3,A,T12,1X,I10))")&
        "info    = ",info
      IF (info > 0) THEN
        WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
          "ifail   = ",ifail
      END IF
      CALL stop_program(routineP,"Error in dsyevx")
    END IF

!   *** Release work storage ***

    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work")

#endif
    eigenvalues(1:neig) = w(1:neig)
    DEALLOCATE (w,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_syevx


! *****************************************************************************

  SUBROUTINE cp_fm_power(matrix,work,exponent,threshold,n_dependent,&
                                work_syevx)

!   Purpose: Raise the real symmetric n by n matrix to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched.

!   History: - Creation (29.03.1999, Matthias Krack)
!            - Parallelised using BLACS and ScaLAPACK (06.06.2001, MK)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                 :: matrix,work
    REAL(KIND = dp), INTENT(IN)                      :: exponent,threshold,work_syevx
    INTEGER, INTENT(OUT)                      :: n_dependent

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE cp_fm_power (MODULE cp_fm_diag)"

!   *** Local variables ***

    REAL(KIND = dp) :: f,p
    INTEGER  :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
                irow_local,istat,mypcol,myprow,ncol_block,ncol_global,npcol,&
                nprow,nrow_block,nrow_global

    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: eigenvalues
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: eigenvectors
    TYPE(cp_blacs_env_type), POINTER    :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_power","I","",handle)


    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    n_dependent = 0
    p = 0.5_dp*exponent

    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global

    ALLOCATE (eigenvalues(ncol_global),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues",ncol_global*dp_size)
    eigenvalues(:) = 0.0_dp

!   *** Compute the eigenvectors and eigenvalues ***

    CALL cp_fm_syevx(matrix,work,eigenvalues,ncol_global,work_syevx)
    IF (matrix%matrix_struct%para_env%mepos.eq.0) THEN
      write(cp_logger_get_default_unit_nr(),'(T2,A,2E16.8)') &
                  "S inversion : smallest and largest evals",eigenvalues(1),eigenvalues(ncol_global)
      IF (eigenvalues(1).LE.0.0_dp) THEN
        write(cp_logger_get_default_unit_nr(),'(T2,A)') &
                    "S inversion : warning : S has a negative eigenvalue, tighten EPS_DEFAULT"
      ENDIF
      write(cp_logger_get_default_unit_nr(),'(T2,A,1E16.8)') &
                    "S inversion : condition number of S = ",eigenvalues(ncol_global)/eigenvalues(1)
    ENDIF

#if defined(__SCALAPACK)

    nrow_block = work%matrix_struct%nrow_block
    ncol_block = work%matrix_struct%ncol_block

    eigenvectors => work%local_data

!   *** Build matrix**exponent with eigenvector quenching ***

    p = 0.5_dp*exponent

    n_dependent = 0

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1

        ipcol = indxg2p(icol_global,ncol_block,mypcol,&
             work%matrix_struct%first_p_pos(2),npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,&
               work%matrix_struct%first_p_pos(2),npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,&
                 work%matrix_struct%first_p_pos(1),nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,&
                   work%matrix_struct%first_p_pos(1),nprow)
              eigenvectors(irow_local,icol_local) = 0.0_dp
            END IF
          END DO
        END IF

      ELSE

        f = eigenvalues(icol_global)**p

        ipcol = indxg2p(icol_global,ncol_block,mypcol,&
             work%matrix_struct%first_p_pos(2),npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,&
               work%matrix_struct%first_p_pos(2),npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,&
                 work%matrix_struct%first_p_pos(1),nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,&
                   work%matrix_struct%first_p_pos(1),nprow)
              eigenvectors(irow_local,icol_local) =&
                f*eigenvectors(irow_local,icol_local)
            END IF
          END DO
        END IF

      END IF

    END DO

#else

    eigenvectors => work%local_data

!   *** Build matrix**exponent with eigenvector quenching ***

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1
        eigenvectors(1:nrow_global,icol_global) = 0.0_dp

      ELSE

        f = eigenvalues(icol_global)**p
        eigenvectors(1:nrow_global,icol_global) =&
          f*eigenvectors(1:nrow_global,icol_global)

      END IF

    END DO

#endif
    CALL cp_fm_syrk("U","N",ncol_global,1.0_dp,work,1,1,0.0_dp,matrix)

    DEALLOCATE (eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_power

!   ***************************************************************************

  SUBROUTINE cp_fm_jacobi(matrix,eigenvectors,eigenvalues,epsi,Anz_max,&
                          thresh,blocksize)

!   Purpose: Calculates Eigenvalues und Eigenvectors from full symmetric matrix

!   History: - Creation (16.09.2002, Martin Fengler) Has its origin in cp_fm_syevx

!   Remark: The Eigenvalues and Eigenvectors are NOT sorted! In our purpose it
!           is not necessary to sort them since they are already at their
!           correct posistions after a couple of SCF-Steps.
!           So this routine should only be used after some general eigenvalue
!           solvers worked on the same problem

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
    REAL(KIND=dp), INTENT(IN)                :: epsi
    INTEGER, INTENT(IN)                      :: Anz_max
    REAL(KIND=dp), INTENT(IN)                :: thresh
    INTEGER, INTENT(IN), OPTIONAL            :: blocksize

    INTEGER                                  :: Anz, bs, handle, I, J, N, p, q
    REAL(KIND=dp)                            :: A_off, c, J_pp, J_pq, J_qp, &
                                                J_qq, s, tan_theta, tau
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: c_ip, c_iq
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: A, EV

! Gibt an wieviele Eigenwerte/vektoren berechnet werden solln
! Gibt die maximal Anzahl an Jacobi-Iterationen an
! epsi: Abbruchkriterium, wenn epsi kleiner als die Quadratsumme
! der Ausserdiagonalwerte ist
! thresh: Ein Element wird nur rotiert, wenn es groesser als thresh ist.
! Enthaelt die Eigenwerte 1:blocksize (vom kleinsten bis zum groessten)
! Da ab einem gewissen SCF-SChritt nicht mehr umsortiert werden muss
! EV: Eigenvektoren
!   ***** LOKALE VARIABLEN *******
!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_jacobi","I","",handle)

    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die diagonisiert werden soll


    ALLOCATE(c_ip(N),c_iq(N))   ! Speicher fuer den lokalen Eigenwertvektor

    eigenvalues(:) = 0.0_dp      ! Contains at the end the eigenvalues (more precise: The diagonal of A)
    A => matrix%local_data       ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell

!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***

       ! Arguments checken
       IF(.NOT. PRESENT(blocksize)) THEN
         bs = N
       ELSE
         bs = blocksize
       ENDIF


       tan_theta = 0.0_dp
       tau       = 0.0_dp

       ! Q zur Einheitsmatrix machen
       EV = 0.0_dp 
       FORALL (I=1:N)
         EV(I,I)=1.0_dp
       END FORALL

       ! Summe der Off-Diagonalwerte bis BLocksize ausrechnen
       A_off = 0.0_dp
       DO J=2,bs
         DO I=1,(J-1)
           A_off = A_off + A(I,J)*A(I,J)
         ENDDO
       ENDDO
       A_off=2.0_dp * A_off

      ! Anzahl der Durchgaenge
       Anz=1

      DO
         IF (A_off < epsi ) THEN
           EXIT
         ENDIF


         DO q=2,bs
          DO p=1,q-1


          IF(ABS(A(p,q))>thresh) THEN

           ! Durch das Wegrotieren erniedirgt sich die Quadratsumme der
           ! Ausserdiagonalwerte:
          A_off = A_off - 2.0_dp*A(p,q)*A(p,q)

          tau = (A(q,q)-A(p,p))/(2.0_dp*A(p,q))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Schnelle Implementierung von: J' * A * J

          ! Erst wird von links dranmultiplitziert
          ! J'*A , weil die Symmetrie mir die Strides rettet

          J_pp=c  
          J_pq=s
          J_qp=-s 
          J_qq=c


          c_ip(1:bs)=A(1:bs,p)*J_pp + A(1:bs,q)*J_qp
          c_iq(1:bs)=A(1:bs,p)*J_pq + A(1:bs,q)*J_qq

          A(1:bs,p)=c_ip(1:bs)
          A(1:bs,q)=c_iq(1:bs)


          ! Jetzt dann von rechts
          ! A*J

          A(p,p) = J_pp*A(p,p) + J_qp*A(q,p)
          A(q,q) = J_pq*A(p,q) + J_qq*A(q,q)

          A(p,q) = 0.0_dp
          A(q,p) = 0.0_dp

          ! Die 4 Elemente A(q,q),A(p,q),A(q,p),A(q,q) werden unnoetig mit transferiert...soll mal latte sein
           A(p,1:bs)=A(1:bs,p)
           A(q,1:bs)=A(1:bs,q)


          !A(p,1:p-1)  = A(1:p-1,p)
          !A(p,p+1:q-1)= A(p+1:q-1,p)
          !A(p,q+1:bs) = A(q+1:bs,p)

          !A(q,1:p-1)  = A(1:p-1,q)
          !A(q,p+1:q-1)= A(p+1:q-1,q)
          !A(q,q+1:bs) = A(q+1:bs,q)



          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Erst wird von rechts dranmultiplitziert
          ! Q * J
          IF (Anz==1) THEN
            c_ip(1:q)=EV(1:q,p)*J_pp + EV(1:q,q)*J_qp
            c_iq(1:q)=EV(1:q,p)*J_pq + EV(1:q,q)*J_qq

            EV(1:q,p)=c_ip
            EV(1:q,q)=c_iq
          ELSE
            c_ip=EV(:,p)*J_pp + EV(:,q)*J_qp
            c_iq=EV(:,p)*J_pq + EV(:,q)*J_qq

            EV(:,p)=c_ip
            EV(:,q)=c_iq
          ENDIF

        ENDIF


        IF(A_off<epsi) THEN
          EXIT
        ENDIF


      ENDDO  ! for-loop Innere

      IF(A_off<epsi) THEN
       EXIT
      ENDIF


    ENDDO   ! for-loop

    Anz=Anz+1

    IF(Anz>Anz_max) THEN
      EXIT
    ENDIF

   ENDDO  ! END-WHILE

  IF(A_off<=epsi) THEN
    WRITE(*,*) 'Die Quadrat-Summe der Ausserdiagnonalelemente ist mit unter epsi: ',A_off
  ENDIF

 ! Abschliessend noch sortieren der Eigenwerte und Eigenvektoren
 ! nach der Groesse

 WRITE(*,*) 'Anzahl der Zyklischen Durchlaeufe: ',(Anz-1)

 FORALL (I=1:blocksize)
    eigenvalues(I)=A(I,I)
 END FORALL


!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***


!   *** Release work storage ***
    DEALLOCATE (c_ip,c_iq)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_jacobi

! *****************************************************************************

  SUBROUTINE cp_fm_block_jacobi_scaled(matrix,eigenvectors,eigenvalues,thresh,start_sec_block)
!
!   Purpose: Calculates Block diagonalizazion from full symmetric matrix
!   History: - Creation (07.10.2002, Martin Fengler) Has its origin in cp_fm_syevx
!   This Routines scales on-the-fly the matrix elements to detect significant influences on the eigenvalues
!   IT DOES ONLY ONE SWEEP!

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: eigenvalues
    REAL(KIND=dp), INTENT(IN)                :: thresh
    INTEGER, INTENT(IN)                      :: start_sec_block

    INTEGER                                  :: handle, I, N, p, q
    REAL(KIND=dp)                            :: c, k, max_thresh, new_thresh, &
                                                s, tan_theta, tau, temp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: c_ip, c_iq
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: A, EV

! Gibt an wo der zweite Block BEGINNT
! thresh: Ein SKALIERTES Element wird nur rotiert, wenn es grsser als thresh ist.
! Enthaelt die Eigenwerte 1:N(vom kleinsten bis zum grssten)
! Da ab einem gewissen SCF-SChritt nicht mehr umsortiert werden muss
! EV: Eigenvektoren
!   ***** LOKALE VARIABLEN *******
!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_block_jacobi_scaled","I","",handle)

    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die diagonisiert werden soll

    ALLOCATE(c_ip(N),c_iq(N))   ! Speicher fuer den lokalen Eigenwertvektor

    A => matrix%local_data       ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell


!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***

       tan_theta = 0.0_dp
       tau       = 0.0_dp

       ! Q zur Einheitsmatrix machen
       EV = 0.0_dp 
       FORALL (I=1:N)
         EV(I,I)=1.0_dp
       END FORALL

    k=thresh
    max_thresh=0.0_dp
    DO q=start_sec_block,N
      DO p=1,(start_sec_block-1)
       temp=ABS(1.0_dp/(eigenvalues(p)-eigenvalues(q))*A(p,q))
       IF(temp>max_thresh) THEN
         max_thresh=temp
       ENDIF
      ENDDO
    ENDDO

    new_thresh=k*max_thresh


    DO q=start_sec_block,N
      DO p=1,(start_sec_block-1)

          temp=ABS(1.0_dp/(eigenvalues(p)-eigenvalues(q))*A(p,q))
          IF(temp>new_thresh .OR. temp>0.00001_dp) THEN

           ! Durch das Wegrotieren erniedirgt sich die Quadratsumme der
           ! Ausserdiagonalwerte:

          tau = (A(q,q)-A(p,p))/(2.0_dp*A(p,q))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
            c_ip = c*EV(:,p) - s*EV(:,q)
            c_iq = s*EV(:,p) + c*EV(:,q)

            EV(:,p)=c_ip
            EV(:,q)=c_iq

         ENDIF


      ENDDO  ! for-loop Innere
    ENDDO   ! for-loop

!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***


!   *** Release work storage ***
    DEALLOCATE (c_ip,c_iq)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_block_jacobi_scaled

! *****************************************************************************


  SUBROUTINE cp_fm_block_jacobi_classic(matrix,eigenvectors,eigval,thresh,&
                                        start_sec_block)
!
!   Purpose: Calculates Block diagonalizazion from full symmetric matrix
!   History: - Creation (07.10.2002, Martin Fengler) Has its origin in cp_fm_syevx
!   This Routines rotates only elements which are larger than a threshold "thresh"
!   IT DOES ONLY ONE SWEEP!

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER        :: eigenvectors,matrix
    REAL(KIND = dp), INTENT(IN)                      :: eigval(:)
    INTEGER, INTENT(IN)                       :: start_sec_block  ! Gibt an wo der zweite Block BEGINNT
    REAL(KIND = dp), INTENT(IN)                      :: thresh ! thresh: Ein Element wird nur rotiert, wenn es groesser als thresh ist.

!   *** Local variables ***
    INTEGER :: handle
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: A,EV         ! EV: Eigenvektoren

!   ***** LOKALE VARIABLEN *******

    REAL(KIND = dp) :: tan_theta,tau,c,s
    INTEGER  :: I,J,q,p,bs,N,q_loc
    REAL(KIND = dp), DIMENSION (:),ALLOCATABLE :: c_ip

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: numroc,indxl2g,indxg2l
    INTEGER :: myprow,mypcol,nprow,npcol,ictxt_loc,block_dim_row,block_dim_col,info,ev_row_block_size,iam,counter,status,ierr
    INTEGER :: source,allgrp,mynumrows,ictxt,temp_row,temp_col,temp_nprow,temp_npcol,mype,npe,a_row,a_col,a_nb,a_mb,iglob,jglob,&
               i_loc,j_loc
    INTEGER, DIMENSION(9) :: desca,descz,desc_a_block, desc_ev_loc
    REAL(KIND = dp), ALLOCATABLE :: A_loc(:,:), EV_loc(:,:),diagonale(:)
    TYPE(cp_blacs_env_type), POINTER    :: context
#endif

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_block_jacobi_classic","I","",handle)

#if defined(__SCALAPACK)

    ! 'Umgebung fuer die Kohn-Sham-Matrix und die Eigenvektoren uebernehmen!'

    context =>  matrix%matrix_struct%context
    source =    matrix%matrix_struct%para_env%source ! globenv%ionode
    allgrp =    matrix%matrix_struct%para_env%group ! globenv%ionode

    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    N = matrix%matrix_struct%nrow_global

    A => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)
    EV => eigenvectors%local_data
    descz(:) = eigenvectors%matrix_struct%descriptor(:)

! Kopiere den Block, der wegrotiert werden soll zunaechst auf den Masterprozessor, und anschliessend
! per Broadcast an alle Prozis
! ACHTUNG start_sec_block sagt aus WO der ZWEITE Block STARTET!!!
! Der Block wird mitsamt dem OO-Block bearbeitet

    block_dim_row = start_sec_block-1
    block_dim_col = N - block_dim_row
    ALLOCATE(A_loc(block_dim_row,block_dim_col))

    mype=matrix%matrix_struct%para_env%mepos 
    npe=matrix%matrix_struct%para_env%num_pe 
    ictxt    =matrix%matrix_struct%context%group
    ! get a new context
    ictxt_loc=matrix%matrix_struct%para_env%group
    CALL cp_blacs_gridinit(ictxt_loc,'Row-major',NPROW*NPCOL,1)

    CALL descinit(desc_a_block,block_dim_row,block_dim_col,block_dim_row,&
                  block_dim_col,0,0,ictxt_loc,block_dim_row,info)

    CALL pdgemr2d(block_dim_row,block_dim_col,A,1,start_sec_block,desca,&
                  A_loc,1,1,desc_a_block,ictxt)
    ! Jetzt sind eigentlich nur im Master-Prozess Daten reingekommen
    CALL mp_bcast(A_loc,0,allgrp)

    ! Da nun jeder ueber den oberen Block verfuegt, koennen wir jetzt die Eigenvektoren so umsortieren, dass
    ! man ein NN*1 Prozessgrid hat, und somit jeder Prozessor ueber ein Buendel von Zeilel verfuegt,
    ! die selbst unabhaengig modifizieren darf.

    ! Aufsetzen der Eigenvektorverteilung
    iam = mype
    ev_row_block_size = n/(nprow*npcol)
    mynumrows = NUMROC(N,ev_row_block_size,iam,0,NPROW*NPCOL)

    ALLOCATE(EV_loc(mynumrows,N),c_ip(mynumrows))

    CALL descinit(desc_ev_loc,N,N,ev_row_block_size,N,0,0,ictxt_loc,&
                  mynumrows,info)

    CALL pdgemr2d(N,N,EV,1,1,descz,EV_loc,1,1,desc_ev_loc,ictxt)

!   *** START Diagonalising matrix ***

    ! Eigentliche Blockdiagonalisierung

    q_loc = 0

    DO q=start_sec_block,N
      q_loc = q_loc + 1
      DO p=1,(start_sec_block-1)

          IF(ABS(A_loc(p,q_loc))>thresh) THEN

          tau = (eigval(q)-eigval(p))/(2.0_dp*A_loc(p,q_loc))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV_loc(:,p) - s*EV_loc(:,q)
          !  c_iq = s*EV_loc(:,p) + c*EV_loc(:,q)

          !  EV(:,p)=c_ip
          !  EV(:,q)=c_iq

          CALL dcopy(mynumrows,EV_loc(1,p),1,c_ip(1),1)
          CALL dscal(mynumrows,c,EV_loc(1,p),1)
          CALL daxpy(mynumrows,-s,EV_loc(1,q),1,EV_loc(1,p),1)
          CALL dscal(mynumrows,c,EV_loc(1,q),1)
          CALL daxpy(mynumrows,s,c_ip(1),1,EV_loc(1,q),1)

         END IF

      END DO
    END DO

!   *** END Diagonalising matrix ***

    ! Nun muessen die Eigenvektoren wieder in die alte Verteilung zurueckverschickt werden.
    ! Verschicke EVs zurueck!'
    CALL pdgemr2d(N,N,EV_loc,1,1,desc_ev_loc,EV,1,1,descz,ictxt)

    ! Speicher freigeben
    DEALLOCATE(A_loc,EV_loc,c_ip)

    CALL cp_blacs_gridexit(ictxt_loc)

    ! ...und tschuess!

#else

    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die bearbeitet werden soll


    ALLOCATE(c_ip(N))   ! Speicher fuer den lokalen Eigenwertvektor

    A => matrix%local_data       ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell

!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***

    tan_theta = 0.0_dp
    tau       = 0.0_dp

    DO q=start_sec_block,N
      DO p=1,(start_sec_block-1)

          IF(ABS(A(p,q))>thresh) THEN

          tau = (eigval(q)-eigval(p))/(2.0_dp*A(p,q))

          tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau))

          ! cos theta
          c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
          s = tan_theta*c

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV(:,p) - s*EV(:,q)
          !  c_iq = s*EV(:,p) + c*EV(:,q)

          !  EV(:,p)=c_ip
          !  EV(:,q)=c_iq

          CALL dcopy(N,EV(1,p),1,c_ip(1),1)
          CALL dscal(N,c,EV(1,p),1)
          CALL daxpy(N,-s,EV(1,q),1,EV(1,p),1)
          CALL dscal(N,c,EV(1,q),1)
          CALL daxpy(N,s,c_ip(1),1,EV(1,q),1)

         ENDIF


      ENDDO  ! for-loop Innere
    ENDDO   ! for-loop

!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***


!   *** Release work storage ***
    DEALLOCATE (c_ip)

#endif

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_block_jacobi_classic

! *****************************************************************************

END MODULE cp_fm_diag
