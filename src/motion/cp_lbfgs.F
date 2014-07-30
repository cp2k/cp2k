! *****************************************************************************
!> \brief LBFGS-B routine (version 3.0, April 25, 2011)
!> \note
!>      L-BFGS-B (version 3.0, April 25, 2011) converted to Fortran 90 module
!> \par History
!>      02.2005 Update to the new version 2.4 and deleting the blas part of
!>              the code (Teodoro Laino)
!>      11.2012 New version 3.0 converted to Fortran 90 (Matthias Krack)
!> \author Fawzi Mohamed (first version)
! *****************************************************************************
MODULE cp_lbfgs
  USE bibliography,                    ONLY: Byrd1995,&
                                             cite_reference
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setulb

CONTAINS

!===========   L-BFGS-B (version 3.0.  April 25, 2011  =================
!
!     This is a modified version of L-BFGS-B. Minor changes in the updat
!     code appear preceded by a line comment as follows
!
!     c-jlm-jn
!
!     Major changes are described in the accompanying paper:
!
!         Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
!         L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constraine
!         Optimization"  (2011). To appear in  ACM Transactions on
!         Mathematical Software,
!
!     The paper describes an improvement and a correction to Algorithm 7
!     It is shown that the performance of the algorithm can be improved
!     significantly by making a relatively simple modication to the subs
!     minimization phase. The correction concerns an error caused by the
!     of routine dpmeps to estimate machine precision.
!
!     The total work space **wa** required by the new version is
!
!                  2*m*n + 11m*m + 5*n + 8*m
!
!     the old version required
!
!                  2*m*n + 12m*m + 4*n + 12*m
!
!
!            J. Nocedal  Department of Electrical Engineering and
!                        Computer Science.
!                        Northwestern University. Evanston, IL. USA
!
!
!           J.L Morales  Departamento de Matematicas,
!                        Instituto Tecnologico Autonomo de Mexico
!                        Mexico D.F. Mexico.
!
!                        March  2011
!
!=======================================================================
! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param x ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param f ...
!> \param g ...
!> \param factr ...
!> \param pgtol ...
!> \param wa ...
!> \param iwa ...
!> \param task ...
!> \param iprint ...
!> \param csave ...
!> \param lsave ...
!> \param isave ...
!> \param dsave ...
! *****************************************************************************
      SUBROUTINE setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,&
                       task, iprint, csave, lsave, isave, dsave)

    INTEGER                                  :: n, m
    REAL(KIND=dp)                            :: x(n), l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp)                            :: f, g(n), factr, pgtol, &
                                                wa(2*m*n + 5*n + 11*m*m + 8*m)
    INTEGER                                  :: iwa(3*n)
    CHARACTER(LEN=60)                        :: task
    INTEGER                                  :: iprint
    CHARACTER(LEN=60)                        :: csave
    LOGICAL                                  :: lsave(4)
    INTEGER                                  :: isave(44)
    REAL(KIND=dp)                            :: dsave(29)

    INTEGER                                  :: ld, lr, lsnd, lss, lsy, lt, &
                                                lwa, lwn, lws, lwt, lwy, lxp, &
                                                lz

!
!-jlm-jn
!     ************
!
!     Subroutine setulb
!
!     This subroutine partitions the working arrays wa and iwa, and
!       then uses the limited memory BFGS method to solve the bound
!       constrained optimization problem by calling mainlb.
!       (The direct method will be used in the subspace minimization.)
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     x is a double precision array of dimension n.
!       On entry x is an approximation to the solution.
!       On exit x is the current approximation.
!
!     l is a double precision array of dimension n.
!       On entry l is the lower bound on x.
!       On exit l is unchanged.
!
!     u is a double precision array of dimension n.
!       On entry u is the upper bound on x.
!       On exit u is unchanged.
!
!     nbd is an integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds, and
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     f is a double precision variable.
!       On first entry f is unspecified.
!       On final exit f is the value of the function at x.
!
!     g is a double precision array of dimension n.
!       On first entry g is unspecified.
!       On final exit g is the value of the gradient at x.
!
!     factr is a double precision variable.
!       On entry factr >= 0 is specified by the user.  The iteration
!         will stop when
!
!         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!
!         where epsmch is the machine precision, which is automatically
!         generated by the code. Typical values for factr: 1.d+12 for
!         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
!         high accuracy.
!       On exit factr is unchanged.
!
!     pgtol is a double precision variable.
!       On entry pgtol >= 0 is specified by the user.  The iteration
!         will stop when
!
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!
!         where pg_i is the ith component of the projected gradient.
!       On exit pgtol is unchanged.
!
!     wa is a double precision working array of length
!       (2mmax + 5)nmax + 12mmax^2 + 12mmax.
!
!     iwa is an integer working array of length 3nmax.
!
!     task is a working string of characters of length 60 indicating
!       the current job when entering and quitting this subroutine.
!
!     iprint is an integer variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     csave is a working string of characters of length 60.
!
!     lsave is a logical working array of dimension 4.
!       On exit with 'task' = NEW_X, the following information is
!                                                             available:
!         If lsave(1) = .true.  then  the initial X has been replaced by
!                                     its projection in the feasible set
!         If lsave(2) = .true.  then  the problem is constrained;
!         If lsave(3) = .true.  then  each variable has upper and lower
!                                     bounds;
!
!     isave is an integer working array of dimension 44.
!       On exit with 'task' = NEW_X, the following information is
!                                                             available:
!         isave(22) = the total number of intervals explored in the
!                         search of Cauchy points;
!         isave(26) = the total number of skipped BFGS updates before
!                         the current iteration;
!         isave(30) = the number of current iteration;
!         isave(31) = the total number of BFGS updates prior the current
!                         iteration;
!         isave(33) = the number of intervals explored in the search of
!                         Cauchy point in the current iteration;
!         isave(34) = the total number of function and gradient
!                         evaluations;
!         isave(36) = the number of function value or gradient
!                                  evaluations in the current iteration;
!         if isave(37) = 0  then the subspace argmin is within the box;
!         if isave(37) = 1  then the subspace argmin is beyond the box;
!         isave(38) = the number of free variables in the current
!                         iteration;
!         isave(39) = the number of active constraints in the current
!                         iteration;
!         n + 1 - isave(40) = the number of variables leaving the set of
!                           active constraints in the current iteration;
!         isave(41) = the number of variables entering the set of active
!                         constraints in the current iteration.
!
!     dsave is a double precision working array of dimension 29.
!       On exit with 'task' = NEW_X, the following information is
!                                                             available:
!         dsave(1) = current 'theta' in the BFGS matrix;
!         dsave(2) = f(x) in the previous iteration;
!         dsave(3) = factr*epsmch;
!         dsave(4) = 2-norm of the line search direction vector;
!         dsave(5) = the machine precision epsmch generated by the code;
!         dsave(7) = the accumulated time spent on searching for
!                                                         Cauchy points;
!         dsave(8) = the accumulated time spent on
!                                                 subspace minimization;
!         dsave(9) = the accumulated time spent on line search;
!         dsave(11) = the slope of the line search function at
!                                  the current point of line search;
!         dsave(12) = the maximum relative step length imposed in
!                                                           line search;
!         dsave(13) = the infinity norm of the projected gradient;
!         dsave(14) = the relative step length in the line search;
!         dsave(15) = the slope of the line search function at
!                                 the starting point of the line search;
!         dsave(16) = the square of the 2-norm of the line search
!                                                      direction vector.
!
!     Subprograms called:
!
!       L-BFGS-B Library ... mainlb.
!
!
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
!       limited memory FORTRAN code for solving bound constrained
!       optimization problems'', Tech. Report, NAM-11, EECS Department,
!       Northwestern University, 1994.
!
!       (Postscript files of these papers are available via anonymous
!        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!-jlm-jn

      IF (task == 'START') THEN
         CALL cite_reference(Byrd1995)
         isave(1)  = m*n
         isave(2)  = m**2
         isave(3)  = 4*m**2
                                            ! ws      m*n
         isave(4)  = 1
                                            ! wy      m*n
         isave(5)  = isave(4)  + isave(1)
                                            ! wsy     m**2
         isave(6)  = isave(5)  + isave(1)
                                            ! wss     m**2
         isave(7)  = isave(6)  + isave(2)
                                            ! wt      m**2
         isave(8)  = isave(7)  + isave(2)
                                            ! wn      4*m**2
         isave(9)  = isave(8)  + isave(2)
                                            ! wsnd    4*m**2
         isave(10) = isave(9)  + isave(3)
                                            ! wz      n
         isave(11) = isave(10) + isave(3)
                                            ! wr      n
         isave(12) = isave(11) + n
                                            ! wd      n
         isave(13) = isave(12) + n
                                            ! wt      n
         isave(14) = isave(13) + n
                                            ! wxp     n
         isave(15) = isave(14) + n
                                            ! wa      8*m
         isave(16) = isave(15) + n
      END IF
      lws  = isave(4)
      lwy  = isave(5)
      lsy  = isave(6)
      lss  = isave(7)
      lwt  = isave(8)
      lwn  = isave(9)
      lsnd = isave(10)
      lz   = isave(11)
      lr   = isave(12)
      ld   = isave(13)
      lt   = isave(14)
      lxp  = isave(15)
      lwa  = isave(16)

      CALL mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol,                        &
        wa(lws),wa(lwy),wa(lsy),wa(lss), wa(lwt),                       &
        wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),wa(lxp),           &
        wa(lwa),                                                        &
        iwa(1),iwa(n+1),iwa(2*n+1),task,iprint,                         &
        csave,lsave,isave(22),dsave)

      RETURN

      END SUBROUTINE setulb

!======================= The end of setulb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param x ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param f ...
!> \param g ...
!> \param factr ...
!> \param pgtol ...
!> \param ws ...
!> \param wy ...
!> \param sy ...
!> \param ss ...
!> \param wt ...
!> \param wn ...
!> \param snd ...
!> \param z ...
!> \param r ...
!> \param d ...
!> \param t ...
!> \param xp ...
!> \param wa ...
!> \param index ...
!> \param iwhere ...
!> \param indx2 ...
!> \param task ...
!> \param iprint ...
!> \param csave ...
!> \param lsave ...
!> \param isave ...
!> \param dsave ...
! *****************************************************************************
      SUBROUTINE mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, &
                        sy, ss, wt, wn, snd, z, r, d, t, xp, wa,        &
                        index, iwhere, indx2, task,                     &
                        iprint, csave, lsave, isave, dsave)
    INTEGER                                  :: n, m
    REAL(KIND=dp)                            :: x(n), l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp) :: f, g(n), factr, pgtol, ws(n, m), wy(n, m), sy(m, m), &
      ss(m, m), wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), z(n), r(n), d(n), &
      t(n), xp(n), wa(8*m)
    INTEGER                                  :: INDEX(n), iwhere(n), indx2(n)
    CHARACTER(LEN=60)                        :: task
    INTEGER                                  :: iprint
    CHARACTER(LEN=60)                        :: csave
    LOGICAL                                  :: lsave(4)
    INTEGER                                  :: isave(23)
    REAL(KIND=dp)                            :: dsave(29)

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    CHARACTER(LEN=3)                         :: word
    INTEGER :: col, head, i, iback, ifun, ileave, info, itail, iter, itfile, &
      iupdat, iword, k, nact, nenter, nfgv, nfree, nintol, nseg, nskip
    LOGICAL                                  :: boxed, cnstnd, prjctd, &
                                                updatd, wrk
    REAL(KIND=dp) :: cachyt, cpu1, cpu2, ddot, ddum, dnorm, dr, dtd, epsmch, &
      fold, gd, gdold, lnscht, rr, sbgnrm, sbtime, stp, stpmx, theta, time, &
      time1, time2, tol, xstep

!-jlm-jn
!     ************
!
!     Subroutine mainlb
!
!     This subroutine solves bound constrained optimization problems by
!       using the compact formula of the limited memory BFGS updates.
!
!     n is an integer variable.
!       On entry n is the number of variables.
!       On exit n is unchanged.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric
!          corrections allowed in the limited memory matrix.
!       On exit m is unchanged.
!
!     x is a double precision array of dimension n.
!       On entry x is an approximation to the solution.
!       On exit x is the current approximation.
!
!     l is a double precision array of dimension n.
!       On entry l is the lower bound of x.
!       On exit l is unchanged.
!
!     u is a double precision array of dimension n.
!       On entry u is the upper bound of x.
!       On exit u is unchanged.
!
!     nbd is an integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds,
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     f is a double precision variable.
!       On first entry f is unspecified.
!       On final exit f is the value of the function at x.
!
!     g is a double precision array of dimension n.
!       On first entry g is unspecified.
!       On final exit g is the value of the gradient at x.
!
!     factr is a double precision variable.
!       On entry factr >= 0 is specified by the user.  The iteration
!         will stop when
!
!         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!
!         where epsmch is the machine precision, which is automatically
!         generated by the code.
!       On exit factr is unchanged.
!
!     pgtol is a double precision variable.
!       On entry pgtol >= 0 is specified by the user.  The iteration
!         will stop when
!
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!
!         where pg_i is the ith component of the projected gradient.
!       On exit pgtol is unchanged.
!
!     ws, wy, sy, and wt are double precision working arrays used to
!       store the following information defining the limited memory
!          BFGS matrix:
!          ws, of dimension n x m, stores S, the matrix of s-vectors;
!          wy, of dimension n x m, stores Y, the matrix of y-vectors;
!          sy, of dimension m x m, stores S'Y;
!          ss, of dimension m x m, stores S'S;
!          yy, of dimension m x m, stores Y'Y;
!          wt, of dimension m x m, stores the Cholesky factorization
!                                  of (theta*S'S+LD^(-1)L'); see eq.
!                                  (2.26) in [3].
!
!     wn is a double precision working array of dimension 2m x 2m
!       used to store the LEL^T factorization of the indefinite matrix
!                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!
!       where     E = [-I  0]
!                     [ 0  I]
!
!     snd is a double precision working array of dimension 2m x 2m
!       used to store the lower triangular part of
!                 N = [Y' ZZ'Y   L_a'+R_z']
!                     [L_a +R_z  S'AA'S   ]
!
!     z(n),r(n),d(n),t(n), xp(n),wa(8*m) are double precision working ar
!       z  is used at different times to store the Cauchy point and
!          the Newton point.
!       xp is used to safeguard the projected Newton direction
!
!     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays.
!
!     index is an integer working array of dimension n.
!       In subroutine freev, index is used to store the free and fixed
!          variables at the Generalized Cauchy Point (GCP).
!
!     iwhere is an integer working array of dimension n used to record
!       the status of the vector x for GCP computation.
!       iwhere(i)=0 or -3 if x(i) is free and has bounds,
!                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
!                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
!                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
!                -1       if x(i) is always free, i.e., no bounds on it.
!
!     indx2 is an integer working array of dimension n.
!       Within subroutine cauchy, indx2 corresponds to the array iorder.
!       In subroutine freev, a list of variables entering and leaving
!       the free set is stored in indx2, and it is passed on to
!       subroutine formk with this information.
!
!     task is a working string of characters of length 60 indicating
!       the current job when entering and leaving this subroutine.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     csave is a working string of characters of length 60.
!
!     lsave is a logical working array of dimension 4.
!
!     isave is an integer working array of dimension 23.
!
!     dsave is a double precision working array of dimension 29.
!
!
!     Subprograms called
!
!       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk,
!
!        errclb, prn1lb, prn2lb, prn3lb, active, projgr,
!
!        freev, cmprlb, matupd, formt.
!
!       Minpack2 Library ... timer
!
!       Linpack Library ... dcopy, ddot.
!
!
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!       Subroutines for Large Scale Bound Constrained Optimization''
!       Tech. Report, NAM-11, EECS Department, Northwestern University,
!       1994.
!
!       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
!       Quasi-Newton Matrices and their use in Limited Memory Methods'',
!       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
!
!       (Postscript files of these papers are available via anonymous
!        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      IF (task == 'START') THEN

         epsmch = EPSILON(one)

         CALL timer(time1)

!        Initialize counters and scalars when task='START'.

!           for the limited memory BFGS matrices:
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .FALSE.
         iback  = 0
         itail  = 0
         iword  = 0
         nact   = 0
         ileave = 0
         nenter = 0
         fold   = zero
         dnorm  = zero
         cpu1   = zero
         gd     = zero
         stpmx  = zero
         sbgnrm = zero
         stp    = zero
         gdold  = zero
         dtd    = zero

!           for operation counts:
         iter   = 0
         nfgv   = 0
         nseg   = 0
         nintol = 0
         nskip  = 0
         nfree  = n
         ifun   = 0
!           for stopping tolerance:
         tol = factr*epsmch

!           for measuring running time:
         cachyt = 0
         sbtime = 0
         lnscht = 0

!           'word' records the status of subspace solutions.
         word = '---'

!           'info' records the termination information.
         info = 0

         itfile = 8
         IF (iprint >= 1) THEN
!                                open a summary file 'iterate.dat'
            OPEN (8, file = 'iterate.dat', status = 'unknown')
         END IF

!        Check the input arguments for errors.

         CALL errclb(n,m,factr,l,u,nbd,task,info,k)
         IF (task(1:5) == 'ERROR') THEN
            CALL prn3lb(n,x,f,task,iprint,info,itfile,                  &
                        iter,nfgv,nintol,nskip,nact,sbgnrm,             &
                        zero,nseg,word,iback,stp,xstep,k,               &
                        cachyt,sbtime,lnscht)
            RETURN
         END IF

         CALL prn1lb(n,m,l,u,x,iprint,itfile,epsmch)

!        Initialize iwhere & project x onto the feasible set.

         CALL active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed)

!        The end of the initialization.

      ELSE
!          restore local variables.

         prjctd = lsave(1)
         cnstnd = lsave(2)
         boxed  = lsave(3)
         updatd = lsave(4)

         nintol = isave(1)
         itfile = isave(3)
         iback  = isave(4)
         nskip  = isave(5)
         head   = isave(6)
         col    = isave(7)
         itail  = isave(8)
         iter   = isave(9)
         iupdat = isave(10)
         nseg   = isave(12)
         nfgv   = isave(13)
         info   = isave(14)
         ifun   = isave(15)
         iword  = isave(16)
         nfree  = isave(17)
         nact   = isave(18)
         ileave = isave(19)
         nenter = isave(20)

         theta  = dsave(1)
         fold   = dsave(2)
         tol    = dsave(3)
         dnorm  = dsave(4)
         epsmch = dsave(5)
         cpu1   = dsave(6)
         cachyt = dsave(7)
         sbtime = dsave(8)
         lnscht = dsave(9)
         time1  = dsave(10)
         gd     = dsave(11)
         stpmx  = dsave(12)
         sbgnrm = dsave(13)
         stp    = dsave(14)
         gdold  = dsave(15)
         dtd    = dsave(16)

!        After returning from the driver go to the point where execution
!        is to resume.

         IF (task(1:5) == 'FG_LN') GOTO 666
         IF (task(1:5) == 'NEW_X') GOTO 777
         IF (task(1:5) == 'FG_ST') GOTO 111
         IF (task(1:4) == 'STOP') THEN
            IF (task(7:9) == 'CPU') THEN
!                                          restore the previous iterate.
               CALL dcopy(n,t,1,x,1)
               CALL dcopy(n,r,1,g,1)
               f = fold
            END IF
            GOTO 999
         END IF
      END IF

!     Compute f0 and g0.

      task = 'FG_START'
!          return to the driver to calculate f and g; reenter at 111.
      GOTO 1000
  111 CONTINUE
      nfgv = 1

!     Compute the infinity norm of the (-) projected gradient.

      CALL projgr(n,l,u,nbd,x,g,sbgnrm)

      IF (iprint >= 1) THEN
         WRITE (*,1002) iter,f,sbgnrm
         WRITE (itfile,1003) iter,nfgv,sbgnrm,f
      END IF
      IF (sbgnrm <= pgtol) THEN
!                                terminate the algorithm.
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         GOTO 999
      END IF

! ----------------- the beginning of the loop --------------------------

  222 CONTINUE
      IF (iprint >= 99) WRITE (*,1001) iter + 1
      iword = -1
!
      IF (.NOT. cnstnd .AND. col > 0) THEN
!                                            skip the search for GCP.
         CALL dcopy(n,x,1,z,1)
         wrk = updatd
         nseg = 0
         GOTO 333
      END IF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Compute the Generalized Cauchy Point (GCP).
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CALL timer(cpu1)
      CALL cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z,                     &
                  m,wy,ws,sy,wt,theta,col,head,                         &
                  wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg,             &
                  iprint, sbgnrm, info, epsmch)
      IF (info /= 0) THEN
!         singular triangular system detected; refresh the lbfgs memory.
         IF(iprint >= 1) WRITE (*, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .FALSE.
         CALL timer(cpu2)
         cachyt = cachyt + cpu2 - cpu1
         GOTO 222
      END IF
      CALL timer(cpu2)
      cachyt = cachyt + cpu2 - cpu1
      nintol = nintol + nseg

!     Count the entering and leaving variables for iter > 0;
!     find the index set of free and active variables at the GCP.

      CALL freev(n,nfree,index,nenter,ileave,indx2,                     &
                 iwhere,wrk,updatd,cnstnd,iprint,iter)
      nact = n - nfree

  333 CONTINUE

!     If there are no free variables or B=theta*I, then
!                                        skip the subspace minimization.

      IF (nfree == 0 .OR. col == 0) GOTO 555

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Subspace minimization.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CALL timer(cpu1)

!     Form  the LEL^T factorization of the indefinite
!       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!       where     E = [-I  0]
!                     [ 0  I]

      IF (wrk) CALL formk(n,nfree,index,nenter,ileave,indx2,iupdat,     &
                       updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
      IF (info /= 0) THEN
!          nonpositive definiteness in Cholesky factorization;
!          refresh the lbfgs memory and restart the iteration.
         IF(iprint >= 1) WRITE (*, 1006)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .FALSE.
         CALL timer(cpu2)
         sbtime = sbtime + cpu2 - cpu1
         GOTO 222
      END IF

!        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
!                                                   from 'cauchy').
      CALL cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,                     &
                 theta,col,head,nfree,cnstnd,info)
      IF (info /= 0) GOTO 444

!-jlm-jn   call the direct method.

      CALL subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy,      &
                 theta, x, g, col, head, iword, wa, wn, iprint, info)
  444 CONTINUE
      IF (info /= 0) THEN
!          singular triangular system detected;
!          refresh the lbfgs memory and restart the iteration.
         IF(iprint >= 1) WRITE (*, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .FALSE.
         CALL timer(cpu2)
         sbtime = sbtime + cpu2 - cpu1
         GOTO 222
      END IF

      CALL timer(cpu2)
      sbtime = sbtime + cpu2 - cpu1
  555 CONTINUE

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Line search and optimality tests.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     Generate the search direction d:=z-x.

      DO 40 i = 1, n
         d(i) = z(i) - x(i)
   40 CONTINUE
      CALL timer(cpu1)
  666 CONTINUE
      CALL lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,      &
                  dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,       &
                  boxed,cnstnd,csave,isave(22),dsave(17))
      IF (info /= 0 .OR. iback >= 20) THEN
!          restore the previous iterate.
         CALL dcopy(n,t,1,x,1)
         CALL dcopy(n,r,1,g,1)
         f = fold
         IF (col == 0) THEN
!             abnormal termination.
            IF (info == 0) THEN
               info = -9
!                restore the actual number of f and g evaluations etc.
               nfgv = nfgv - 1
               ifun = ifun - 1
               iback = iback - 1
            END IF
            task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
            iter = iter + 1
            GOTO 999
         ELSE
!             refresh the lbfgs memory and restart the iteration.
            IF(iprint >= 1) WRITE (*, 1008)
            IF (info == 0) nfgv = nfgv - 1
            info   = 0
            col    = 0
            head   = 1
            theta  = one
            iupdat = 0
            updatd = .FALSE.
            task   = 'RESTART_FROM_LNSRCH'
            CALL timer(cpu2)
            lnscht = lnscht + cpu2 - cpu1
            GOTO 222
         END IF
      ELSE IF (task(1:5) == 'FG_LN') THEN
!          return to the driver for calculating f and g; reenter at 666.
         GOTO 1000
      ELSE
!          calculate and print out the quantities related to the new X.
         CALL timer(cpu2)
         lnscht = lnscht + cpu2 - cpu1
         iter = iter + 1

!        Compute the infinity norm of the projected (-)gradient.

         CALL projgr(n,l,u,nbd,x,g,sbgnrm)

!        Print iteration information.

         CALL prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,              &
                     sbgnrm,nseg,word,iword,iback,stp,xstep)
         GOTO 1000
      END IF
  777 CONTINUE

!     Test for termination.

      IF (sbgnrm <= pgtol) THEN
!                                terminate the algorithm.
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         GOTO 999
      END IF

      ddum = MAX(ABS(fold), ABS(f), one)
      IF ((fold - f) <= tol*ddum) THEN
!                                        terminate the algorithm.
         task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
         IF (iback >= 10) info = -5
!           i.e., to issue a warning if iback>10 in the line search.
         GOTO 999
      END IF

!     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.

      DO 42 i = 1, n
         r(i) = g(i) - r(i)
   42 CONTINUE
      rr = ddot(n,r,1,r,1)
      IF (stp == one) THEN
         dr = gd - gdold
         ddum = -gdold
      ELSE
         dr = (gd - gdold)*stp
         CALL dscal(n,stp,d,1)
         ddum = -gdold*stp
      END IF

      IF (dr <= epsmch*ddum) THEN
!                            skip the L-BFGS update.
         nskip = nskip + 1
         updatd = .FALSE.
         IF (iprint >= 1) WRITE (*,1004) dr, ddum
         GOTO 888
      END IF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Update the L-BFGS matrix.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      updatd = .TRUE.
      iupdat = iupdat + 1

!     Update matrices WS and WY and form the middle matrix in B.

      CALL matupd(n,m,ws,wy,sy,ss,d,r,itail,                            &
                  iupdat,col,head,theta,rr,dr,stp,dtd)

!     Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
!        Store T in the upper triangular of the array wt;
!        Cholesky factorize T to J*J' with
!           J' stored in the upper triangular of wt.

      CALL formt(m,wt,sy,ss,col,theta,info)

      IF (info /= 0) THEN
!          nonpositive definiteness in Cholesky factorization;
!          refresh the lbfgs memory and restart the iteration.
         IF(iprint >= 1) WRITE (*, 1007)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .FALSE.
         GOTO 222
      END IF

!     Now the inverse of the middle matrix in B is

!       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
!       [ -L*D^(-1/2)   J ] [  0        J'          ]

  888 CONTINUE

! -------------------- the end of the loop -----------------------------

      GOTO 222
  999 CONTINUE
      CALL timer(time2)
      time = time2 - time1
      CALL prn3lb(n,x,f,task,iprint,info,itfile,                        &
                  iter,nfgv,nintol,nskip,nact,sbgnrm,                   &
                  time,nseg,word,iback,stp,xstep,k,                     &
                  cachyt,sbtime,lnscht)
 1000 CONTINUE

!     Save local variables.

      lsave(1)  = prjctd
      lsave(2)  = cnstnd
      lsave(3)  = boxed
      lsave(4)  = updatd

      isave(1)  = nintol
      isave(3)  = itfile
      isave(4)  = iback
      isave(5)  = nskip
      isave(6)  = head
      isave(7)  = col
      isave(8)  = itail
      isave(9)  = iter
      isave(10) = iupdat
      isave(12) = nseg
      isave(13) = nfgv
      isave(14) = info
      isave(15) = ifun
      isave(16) = iword
      isave(17) = nfree
      isave(18) = nact
      isave(19) = ileave
      isave(20) = nenter

      dsave(1)  = theta
      dsave(2)  = fold
      dsave(3)  = tol
      dsave(4)  = dnorm
      dsave(5)  = epsmch
      dsave(6)  = cpu1
      dsave(7)  = cachyt
      dsave(8)  = sbtime
      dsave(9)  = lnscht
      dsave(10) = time1
      dsave(11) = gd
      dsave(12) = stpmx
      dsave(13) = sbgnrm
      dsave(14) = stp
      dsave(15) = gdold
      dsave(16) = dtd

 1001 FORMAT (//,'ITERATION ',i5)
 1002 FORMAT                                                            &
     &  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 1003 FORMAT (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x,    &
              1p,2(1x,d10.3))
 1004 FORMAT ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' BFGS update SKIPPED')
 1005 FORMAT (/,                                                        &
     &' Singular triangular system detected;',/,                        &
     &'   refresh the lbfgs memory and restart the iteration.')
 1006 FORMAT (/,                                                        &
     &' Nonpositive definiteness in Cholesky factorization in formk;',/,&
     &'   refresh the lbfgs memory and restart the iteration.')
 1007 FORMAT (/,                                                        &
     &' Nonpositive definiteness in Cholesky factorization in formt;',/,&
     &'   refresh the lbfgs memory and restart the iteration.')
 1008 FORMAT (/,                                                        &
     &' Bad direction in the line search;',/,                           &
     &'   refresh the lbfgs memory and restart the iteration.')

      RETURN

      END SUBROUTINE mainlb

!======================= The end of mainlb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param x ...
!> \param iwhere ...
!> \param iprint ...
!> \param prjctd ...
!> \param cnstnd ...
!> \param boxed ...
! *****************************************************************************
      SUBROUTINE active(n, l, u, nbd, x, iwhere, iprint,                &
                        prjctd, cnstnd, boxed)

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp)                            :: x(n)
    INTEGER                                  :: iwhere(n), iprint
    LOGICAL                                  :: prjctd, cnstnd, boxed

    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    INTEGER                                  :: i, nbdd

!     ************
!
!     Subroutine active
!
!     This subroutine initializes iwhere and projects the initial x to
!       the feasible set if necessary.
!
!     iwhere is an integer array of dimension n.
!       On entry iwhere is unspecified.
!       On exit iwhere(i)=-1  if x(i) has no bounds
!                         3   if l(i)=u(i)
!                         0   otherwise.
!       In cauchy, iwhere is given finer gradations.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!     Initialize nbdd, prjctd, cnstnd and boxed.

      nbdd = 0
      prjctd = .FALSE.
      cnstnd = .FALSE.
      boxed = .TRUE.

!     Project the initial x to the easible set if necessary.

      DO 10 i = 1, n
         IF (nbd(i) > 0) THEN
            IF (nbd(i) <= 2 .AND. x(i) <= l(i)) THEN
               IF (x(i) < l(i)) THEN
                  prjctd = .TRUE.
                  x(i) = l(i)
               END IF
               nbdd = nbdd + 1
            ELSE IF (nbd(i) >= 2 .AND. x(i) >= u(i)) THEN
               IF (x(i) > u(i)) THEN
                  prjctd = .TRUE.
                  x(i) = u(i)
               END IF
               nbdd = nbdd + 1
            END IF
         END IF
   10 CONTINUE

!     Initialize iwhere and assign values to cnstnd and boxed.

      DO 20 i = 1, n
         IF (nbd(i) /= 2) boxed = .FALSE.
         IF (nbd(i) == 0) THEN
!                                this variable is always free
            iwhere(i) = -1

!           otherwise set x(i)=mid(x(i), u(i), l(i)).
         ELSE
            cnstnd = .TRUE.
            IF (nbd(i) == 2 .AND. u(i) - l(i) <= zero) THEN
!                   this variable is always fixed
               iwhere(i) = 3
            ELSE
               iwhere(i) = 0
            END IF
         END IF
   20 CONTINUE

      IF (iprint >= 0) THEN
         IF (prjctd) WRITE (*,*)                                        &
     &   'The initial X is infeasible.  Restart with its projection.'
         IF (.NOT. cnstnd)                                              &
            WRITE (*,*) 'This problem is unconstrained.'
      END IF

      IF (iprint > 0) WRITE (*,1001) nbdd

 1001 FORMAT (/,'At X0 ',i9,' variables are exactly at the bounds')

      RETURN

      END SUBROUTINE active

!======================= The end of active =============================

! *****************************************************************************
!> \brief ...
!> \param m ...
!> \param sy ...
!> \param wt ...
!> \param col ...
!> \param v ...
!> \param p ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE bmv(m, sy, wt, col, v, p, info)

    INTEGER                                  :: m
    REAL(KIND=dp)                            :: sy(m, m), wt(m, m)
    INTEGER                                  :: col
    REAL(KIND=dp)                            :: v(2*col), p(2*col)
    INTEGER                                  :: info

    INTEGER                                  :: i, i2, k
    REAL(KIND=dp)                            :: sum

!     ************
!
!     Subroutine bmv
!
!     This subroutine computes the product of the 2m x 2m middle matrix
!       in the compact L-BFGS formula of B and a 2m vector v;
!       it returns the product in p.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     sy is a double precision array of dimension m x m.
!       On entry sy specifies the matrix S'Y.
!       On exit sy is unchanged.
!
!     wt is a double precision array of dimension m x m.
!       On entry wt specifies the upper triangular matrix J' which is
!         the Cholesky factor of (thetaS'S+LD^(-1)L').
!       On exit wt is unchanged.
!
!     col is an integer variable.
!       On entry col specifies the number of s-vectors (or y-vectors)
!         stored in the compact L-BFGS formula.
!       On exit col is unchanged.
!
!     v is a double precision array of dimension 2col.
!       On entry v specifies vector v.
!       On exit v is unchanged.
!
!     p is a double precision array of dimension 2col.
!       On entry p is unspecified.
!       On exit p is the product Mv.
!
!     info is an integer variable.
!       On entry info is unspecified.
!       On exit info = 0       for normal return,
!                    = nonzero for abnormal return when the system
!                                to be solved by dtrsl is singular.
!
!     Subprograms called:
!
!       Linpack ... dtrsl.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      IF (col == 0) RETURN

!     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
!                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

!       solve Jp2=v2+LD^(-1)v1.
      p(col + 1) = v(col + 1)
      DO 20 i = 2, col
         i2 = col + i
         sum = 0.0_dp
         DO 10 k = 1, i - 1
            sum = sum + sy(i,k)*v(k)/sy(k,k)
   10    CONTINUE
         p(i2) = v(i2) + sum
   20 CONTINUE
!     Solve the triangular system
      CALL dtrsl(wt,m,col,p(col+1),11,info)
      IF (info /= 0) RETURN

!       solve D^(1/2)p1=v1.
      DO 30 i = 1, col
         p(i) = v(i)/SQRT(sy(i,i))
   30 CONTINUE

!     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
!                    [  0         J'           ] [ p2 ]   [ p2 ].

!       solve J^Tp2=p2.
      CALL dtrsl(wt,m,col,p(col+1),01,info)
      IF (info /= 0) RETURN

!       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
!                 =-D^(-1/2)p1+D^(-1)L'p2.
      DO 40 i = 1, col
         p(i) = -p(i)/SQRT(sy(i,i))
   40 CONTINUE
      DO 60 i = 1, col
         sum = 0._dp
         DO 50 k = i + 1, col
            sum = sum + sy(k,i)*p(col+k)/sy(i,i)
   50    CONTINUE
         p(i) = p(i) + sum
   60 CONTINUE

      RETURN

      END SUBROUTINE bmv

!======================== The end of bmv ===============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param g ...
!> \param iorder ...
!> \param iwhere ...
!> \param t ...
!> \param d ...
!> \param xcp ...
!> \param m ...
!> \param wy ...
!> \param ws ...
!> \param sy ...
!> \param wt ...
!> \param theta ...
!> \param col ...
!> \param head ...
!> \param p ...
!> \param c ...
!> \param wbp ...
!> \param v ...
!> \param nseg ...
!> \param iprint ...
!> \param sbgnrm ...
!> \param info ...
!> \param epsmch ...
! *****************************************************************************
      SUBROUTINE cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp,  &
                        m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                        v, nseg, iprint, sbgnrm, info, epsmch)
    INTEGER                                  :: n
    REAL(KIND=dp)                            :: x(n), l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp)                            :: g(n)
    INTEGER                                  :: iorder(n), iwhere(n)
    REAL(KIND=dp)                            :: t(n), d(n), xcp(n)
    INTEGER                                  :: m
    REAL(KIND=dp)                            :: sy(m, m), wt(m, m), theta
    INTEGER                                  :: col
    REAL(KIND=dp)                            :: ws(n, col), wy(n, col)
    INTEGER                                  :: head
    REAL(KIND=dp)                            :: p(2*m), c(2*m), wbp(2*m), &
                                                v(2*m)
    INTEGER                                  :: nseg, iprint
    REAL(KIND=dp)                            :: sbgnrm
    INTEGER                                  :: info
    REAL(KIND=dp)                            :: epsmch

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: col2, i, ibkmin, ibp, iter, &
                                                j, nbreak, nfree, nleft, &
                                                pointr
    LOGICAL                                  :: bnded, xlower, xupper
    REAL(KIND=dp)                            :: bkmin, ddot, dibp, dibp2, dt, &
                                                dtm, f1, f2, f2_org, neggi, &
                                                tj, tj0, tl, tsum, tu, wmc, &
                                                wmp, wmw, zibp

!     ************
!
!     Subroutine cauchy
!
!     For given x, l, u, g (with sbgnrm > 0), and a limited memory
!       BFGS matrix B defined in terms of matrices WY, WS, WT, and
!       scalars head, col, and theta, this subroutine computes the
!       generalized Cauchy point (GCP), defined as the first local
!       minimizer of the quadratic
!
!                  Q(x + s) = g's + 1/2 s'Bs
!
!       along the projected gradient direction P(x-tg,l,u).
!       The routine returns the GCP in xcp.
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     x is a double precision array of dimension n.
!       On entry x is the starting point for the GCP computation.
!       On exit x is unchanged.
!
!     l is a double precision array of dimension n.
!       On entry l is the lower bound of x.
!       On exit l is unchanged.
!
!     u is a double precision array of dimension n.
!       On entry u is the upper bound of x.
!       On exit u is unchanged.
!
!     nbd is an integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds, and
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     g is a double precision array of dimension n.
!       On entry g is the gradient of f(x).  g must be a nonzero vector.
!       On exit g is unchanged.
!
!     iorder is an integer working array of dimension n.
!       iorder will be used to store the breakpoints in the piecewise
!       linear path and free variables encountered. On exit,
!         iorder(1),...,iorder(nleft) are indices of breakpoints
!                                which have not been encountered;
!         iorder(nleft+1),...,iorder(nbreak) are indices of
!                                     encountered breakpoints; and
!         iorder(nfree),...,iorder(n) are indices of variables which
!                 have no bound constraits along the search direction.
!
!     iwhere is an integer array of dimension n.
!       On entry iwhere indicates only the permanently fixed (iwhere=3)
!       or free (iwhere= -1) components of x.
!       On exit iwhere records the status of the current x variables.
!       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
!                 0   if x(i) is free and has bounds, and is moved
!                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
!                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
!                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
!                 -1  if x(i) is always free, i.e., it has no bounds.
!
!     t is a double precision working array of dimension n.
!       t will be used to store the break points.
!
!     d is a double precision array of dimension n used to store
!       the Cauchy direction P(x-tg)-x.
!
!     xcp is a double precision array of dimension n used to return the
!       GCP on exit.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     ws, wy, sy, and wt are double precision arrays.
!       On entry they store information that defines the
!                             limited memory BFGS matrix:
!         ws(n,m) stores S, a set of s-vectors;
!         wy(n,m) stores Y, a set of y-vectors;
!         sy(m,m) stores S'Y;
!         wt(m,m) stores the
!                 Cholesky factorization of (theta*S'S+LD^(-1)L').
!       On exit these arrays are unchanged.
!
!     theta is a double precision variable.
!       On entry theta is the scaling factor specifying B_0 = theta I.
!       On exit theta is unchanged.
!
!     col is an integer variable.
!       On entry col is the actual number of variable metric
!         corrections stored so far.
!       On exit col is unchanged.
!
!     head is an integer variable.
!       On entry head is the location of the first s-vector (or y-vector
!         in S (or Y).
!       On exit col is unchanged.
!
!     p is a double precision working array of dimension 2m.
!       p will be used to store the vector p = W^(T)d.
!
!     c is a double precision working array of dimension 2m.
!       c will be used to store the vector c = W^(T)(xcp-x).
!
!     wbp is a double precision working array of dimension 2m.
!       wbp will be used to store the row of W corresponding
!         to a breakpoint.
!
!     v is a double precision working array of dimension 2m.
!
!     nseg is an integer variable.
!       On exit nseg records the number of quadratic segments explored
!         in searching for the GCP.
!
!     sg and yg are double precision arrays of dimension m.
!       On entry sg  and yg store S'g and Y'g correspondingly.
!       On exit they are unchanged.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     sbgnrm is a double precision variable.
!       On entry sbgnrm is the norm of the projected gradient at x.
!       On exit sbgnrm is unchanged.
!
!     info is an integer variable.
!       On entry info is 0.
!       On exit info = 0       for normal return,
!                    = nonzero for abnormal return when the the system
!                              used in routine bmv is singular.
!
!     Subprograms called:
!
!       L-BFGS-B Library ... hpsolb, bmv.
!
!       Linpack ... dscal dcopy, daxpy.
!
!
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!       Subroutines for Large Scale Bound Constrained Optimization''
!       Tech. Report, NAM-11, EECS Department, Northwestern University,
!       1994.
!
!       (Postscript files of these papers are available via anonymous
!        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!     Check the status of the variables, reset iwhere(i) if necessary;
!       compute the Cauchy direction d and the breakpoints t; initialize
!       the derivative f1 and the vector p = W'd (for theta = 1).

      IF (sbgnrm <= zero) THEN
         IF (iprint >= 0) WRITE (*,*) 'Subgnorm = 0.  GCP = X.'
         CALL dcopy(n,x,1,xcp,1)
         RETURN
      END IF
      bnded = .TRUE.
      nfree = n + 1
      nbreak = 0
      ibkmin = 0
      bkmin = zero
      col2 = 2*col
      f1 = zero
      IF (iprint >= 99) WRITE (*,3010)

!     We set p to zero and build it up as we determine d.

      DO 20 i = 1, col2
         p(i) = zero
   20 CONTINUE

!     In the following loop we determine for each variable its bound
!        status and its breakpoint, and update p accordingly.
!        Smallest breakpoint is identified.

      DO 50 i = 1, n
         neggi = -g(i)
         IF (iwhere(i) /= 3 .AND. iwhere(i) /= -1) THEN
!             if x(i) is not a constant and has bounds,
!             compute the difference between x(i) and its bounds.
            IF (nbd(i) <= 2) tl = x(i) - l(i)
            IF (nbd(i) >= 2) tu = u(i) - x(i)

!           If a variable is close enough to a bound
!             we treat it as at bound.
            xlower = nbd(i) <= 2 .AND. tl <= zero
            xupper = nbd(i) >= 2 .AND. tu <= zero

!              reset iwhere(i).
            iwhere(i) = 0
            IF (xlower) THEN
               IF (neggi <= zero) iwhere(i) = 1
            ELSE IF (xupper) THEN
               IF (neggi >= zero) iwhere(i) = 2
            ELSE
               IF (ABS(neggi) <= zero) iwhere(i) = -3
            END IF
         END IF
         pointr = head
         IF (iwhere(i) /= 0 .AND. iwhere(i) /= -1) THEN
            d(i) = zero
         ELSE
            d(i) = neggi
            f1 = f1 - neggi*neggi
!             calculate p := p - W'e_i* (g_i).
            DO 40 j = 1, col
               p(j) = p(j) +  wy(i,pointr)* neggi
               p(col + j) = p(col + j) + ws(i,pointr)*neggi
               pointr = MOD(pointr,m) + 1
   40       CONTINUE
            IF (nbd(i) <= 2 .AND. nbd(i) /= 0                       &
     &                        .AND. neggi < zero) THEN
!                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tl/(-neggi)
               IF (nbreak == 1 .OR. t(nbreak) < bkmin) THEN
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               END IF
            ELSE IF (nbd(i) >= 2 .AND. neggi > zero) THEN
!                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tu/neggi
               IF (nbreak == 1 .OR. t(nbreak) < bkmin) THEN
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               END IF
            ELSE
!                x(i) + d(i) is not bounded.
               nfree = nfree - 1
               iorder(nfree) = i
               IF (ABS(neggi) > zero) bnded = .FALSE.
            END IF
         END IF
   50 CONTINUE

!     The indices of the nonzero components of d are now stored
!       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
!       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.

      IF (theta /= one) THEN
!                   complete the initialization of p for theta not= one.
         CALL dscal(col,theta,p(col+1),1)
      END IF

!     Initialize GCP xcp = x.

      CALL dcopy(n,x,1,xcp,1)

      IF (nbreak == 0 .AND. nfree == n + 1) THEN
!                  is a zero vector, return with the initial xcp as GCP.
         IF (iprint > 100) WRITE (*,1010) (xcp(i), i = 1, n)
         RETURN
      END IF

!     Initialize c = W'(xcp - x) = 0.

      DO 60 j = 1, col2
         c(j) = zero
   60 CONTINUE

!     Initialize derivative f2.

      f2 =  -theta*f1
      f2_org  =  f2
      IF (col > 0) THEN
         CALL bmv(m,sy,wt,col,p,v,info)
         IF (info /= 0) RETURN
         f2 = f2 - ddot(col2,v,1,p,1)
      END IF
      dtm = -f1/f2
      tsum = zero
      nseg = 1
      IF (iprint >= 99)                                               &
         WRITE (*,*) 'There are ',nbreak,'  breakpoints '

!     If there are no breakpoints, locate the GCP and return.

      IF (nbreak == 0) GOTO 888

      nleft = nbreak
      iter = 1


      tj = zero

!------------------- the beginning of the loop -------------------------

  777 CONTINUE

!     Find the next smallest breakpoint;
!       compute dt = t(nleft) - t(nleft + 1).

      tj0 = tj
      IF (iter == 1) THEN
!         Since we already have the smallest breakpoint we need not do
!         heapsort yet. Often only one breakpoint is used and the
!         cost of heapsort is avoided.
         tj = bkmin
         ibp = iorder(ibkmin)
      ELSE
         IF (iter == 2) THEN
!             Replace the already used smallest breakpoint with the
!             breakpoint numbered nbreak > nlast, before heapsort call.
            IF (ibkmin /= nbreak) THEN
               t(ibkmin) = t(nbreak)
               iorder(ibkmin) = iorder(nbreak)
            END IF
!        Update heap structure of breakpoints
!           (if iter=2, initialize heap).
         END IF
         CALL hpsolb(nleft,t,iorder,iter-2)
         tj = t(nleft)
         ibp = iorder(nleft)
      END IF

      dt = tj - tj0

      IF (dt /= zero .AND. iprint >= 100) THEN
         WRITE (*,4011) nseg,f1,f2
         WRITE (*,5010) dt
         WRITE (*,6010) dtm
      END IF

!     If a minimizer is within this interval, locate the GCP and return.

      IF (dtm < dt) GOTO 888

!     Otherwise fix one variable and
!       reset the corresponding component of d to zero.

      tsum = tsum + dt
      nleft = nleft - 1
      iter = iter + 1
      dibp = d(ibp)
      d(ibp) = zero
      IF (dibp > zero) THEN
         zibp = u(ibp) - x(ibp)
         xcp(ibp) = u(ibp)
         iwhere(ibp) = 2
      ELSE
         zibp = l(ibp) - x(ibp)
         xcp(ibp) = l(ibp)
         iwhere(ibp) = 1
      END IF
      IF (iprint >= 100) WRITE (*,*) 'Variable  ',ibp,'  is fixed.'
      IF (nleft == 0 .AND. nbreak == n) THEN
!                                             all n variables are fixed,
!                                                return with xcp as GCP.
         dtm = dt
         GOTO 999
      END IF

!     Update the derivative information.

      nseg = nseg + 1
      dibp2 = dibp**2

!     Update f1 and f2.

!        temporarily set f1 and f2 for col=0.
      f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
      f2 = f2 - theta*dibp2

      IF (col > 0) THEN
!                          update c = c + dt*p.
         CALL daxpy(col2,dt,p,1,c,1)

!           choose wbp,
!           the row of W corresponding to the breakpoint encountered.
         pointr = head
         DO 70 j = 1,col
            wbp(j) = wy(ibp,pointr)
            wbp(col + j) = theta*ws(ibp,pointr)
            pointr = MOD(pointr,m) + 1
   70    CONTINUE

!           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
         CALL bmv(m,sy,wt,col,wbp,v,info)
         IF (info /= 0) RETURN
         wmc = ddot(col2,c,1,v,1)
         wmp = ddot(col2,p,1,v,1)
         wmw = ddot(col2,wbp,1,v,1)

!           update p = p - dibp*wbp.
         CALL daxpy(col2,-dibp,wbp,1,p,1)

!           complete updating f1 and f2 while col > 0.
         f1 = f1 + dibp*wmc
         f2 = f2 + 2.0_dp*dibp*wmp - dibp2*wmw
      END IF

      f2 = MAX(epsmch*f2_org,f2)
      IF (nleft > 0) THEN
         dtm = -f1/f2
         GOTO 777
!                 to repeat the loop for unsearched intervals.
      ELSE IF(bnded) THEN
         f1 = zero
         f2 = zero
         dtm = zero
      ELSE
         dtm = -f1/f2
      END IF

!------------------- the end of the loop -------------------------------

  888 CONTINUE
      IF (iprint >= 99) THEN
         WRITE (*,*)
         WRITE (*,*) 'GCP found in this segment'
         WRITE (*,4010) nseg,f1,f2
         WRITE (*,6010) dtm
      END IF
      IF (dtm <= zero) dtm = zero
      tsum = tsum + dtm

!     Move free variables (i.e., the ones w/o breakpoints) and
!       the variables whose breakpoints haven't been reached.

      CALL daxpy(n,tsum,d,1,xcp,1)

  999 CONTINUE

!     Update c = c + dtm*p = W'(x^c - x)
!       which will be used in computing r = Z'(B(x^c - x) + g).

      IF (col > 0) CALL daxpy(col2,dtm,p,1,c,1)
      IF (iprint > 100) WRITE (*,1010) (xcp(i),i = 1,n)
      IF (iprint >= 99) WRITE (*,2010)

 1010 FORMAT ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
 2010 FORMAT (/,'---------------- exit CAUCHY----------------------',/)
 3010 FORMAT (/,'---------------- CAUCHY entered-------------------')
 4010 FORMAT ('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
 4011 FORMAT (/,'Piece    ',i3,' --f1, f2 at start point ',             &
              1p,2(1x,d11.4))
 5010 FORMAT ('Distance to the next break point =  ',1p,d11.4)
 6010 FORMAT ('Distance to the stationary point =  ',1p,d11.4)

      RETURN

      END SUBROUTINE cauchy

!====================== The end of cauchy ==============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param x ...
!> \param g ...
!> \param ws ...
!> \param wy ...
!> \param sy ...
!> \param wt ...
!> \param z ...
!> \param r ...
!> \param wa ...
!> \param index ...
!> \param theta ...
!> \param col ...
!> \param head ...
!> \param nfree ...
!> \param cnstnd ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index,    &
                       theta, col, head, nfree, cnstnd, info)

    INTEGER                                  :: n, m
    REAL(KIND=dp)                            :: x(n), g(n), ws(n, m), &
                                                wy(n, m), sy(m, m), wt(m, m), &
                                                z(n), r(n), wa(4*m)
    INTEGER                                  :: INDEX(n)
    REAL(KIND=dp)                            :: theta
    INTEGER                                  :: col, head, nfree
    LOGICAL                                  :: cnstnd
    INTEGER                                  :: info

    INTEGER                                  :: i, j, k, pointr
    REAL(KIND=dp)                            :: a1, a2

!     ************
!
!     Subroutine cmprlb
!
!       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
!         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
!
!     Subprograms called:
!
!       L-BFGS-B Library ... bmv.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      IF (.NOT. cnstnd .AND. col > 0) THEN
         DO 26 i = 1, n
            r(i) = -g(i)
   26    CONTINUE
      ELSE
         DO 30 i = 1, nfree
            k = INDEX(i)
            r(i) = -theta*(z(k) - x(k)) - g(k)
   30    CONTINUE
         CALL bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
         IF (info /= 0) THEN
            info = -8
            RETURN
         END IF
         pointr = head
         DO 34 j = 1, col
            a1 = wa(j)
            a2 = theta*wa(col + j)
            DO 32 i = 1, nfree
               k = INDEX(i)
               r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
   32       CONTINUE
            pointr = MOD(pointr,m) + 1
   34    CONTINUE
      END IF

      RETURN

      END SUBROUTINE cmprlb

!======================= The end of cmprlb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param factr ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param task ...
!> \param info ...
!> \param k ...
! *****************************************************************************
      SUBROUTINE errclb(n, m, factr, l, u, nbd, task, info, k)

    INTEGER                                  :: n, m
    REAL(KIND=dp)                            :: factr, l(n), u(n)
    INTEGER                                  :: nbd(n)
    CHARACTER(LEN=60)                        :: task
    INTEGER                                  :: info, k

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i

!     ************
!
!     Subroutine errclb
!
!     This subroutine checks the validity of the input data.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!     Check the input arguments for errors.

      IF (n <= 0) task = 'ERROR: N <= 0'
      IF (m <= 0) task = 'ERROR: M <= 0'
      IF (factr < zero) task = 'ERROR: FACTR < 0'

!     Check the validity of the arrays nbd(i), u(i), and l(i).

      DO 10 i = 1, n
         IF (nbd(i) < 0 .OR. nbd(i) > 3) THEN
!                                                   return
            task = 'ERROR: INVALID NBD'
            info = -6
            k = i
         END IF
         IF (nbd(i) == 2) THEN
            IF (l(i) > u(i)) THEN
!                                    return
               task = 'ERROR: NO FEASIBLE SOLUTION'
               info = -7
               k = i
            END IF
         END IF
   10 CONTINUE

      RETURN

      END SUBROUTINE errclb

!======================= The end of errclb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param nsub ...
!> \param ind ...
!> \param nenter ...
!> \param ileave ...
!> \param indx2 ...
!> \param iupdat ...
!> \param updatd ...
!> \param wn ...
!> \param wn1 ...
!> \param m ...
!> \param ws ...
!> \param wy ...
!> \param sy ...
!> \param theta ...
!> \param col ...
!> \param head ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE formk(n, nsub, ind, nenter, ileave, indx2, iupdat,     &
                       updatd, wn, wn1, m, ws, wy, sy, theta, col,      &
                       head, info)

    INTEGER                                  :: n, nsub, ind(n), nenter, &
                                                ileave, indx2(n), iupdat
    LOGICAL                                  :: updatd
    INTEGER                                  :: m
    REAL(KIND=dp)                            :: wn1(2*m, 2*m), wn(2*m, 2*m), &
                                                ws(n, m), wy(n, m), sy(m, m), &
                                                theta
    INTEGER                                  :: col, head, info

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: col2, dbegin, dend, i, ipntr, &
                                                is, is1, iy, jpntr, js, js1, &
                                                jy, k, k1, m2, pbegin, pend, &
                                                upcl
    REAL(KIND=dp)                            :: ddot, temp1, temp2, temp3, &
                                                temp4

!     ************
!
!     Subroutine formk
!
!     This subroutine forms  the LEL^T factorization of the indefinite
!
!       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!                                                    where E = [-I  0]
!                                                              [ 0  I]
!     The matrix K can be shown to be equal to the matrix M^[-1]N
!       occurring in section 5.1 of [1], as well as to the matrix
!       Mbar^[-1] Nbar in section 5.3.
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     nsub is an integer variable
!       On entry nsub is the number of subspace variables in free set.
!       On exit nsub is not changed.
!
!     ind is an integer array of dimension nsub.
!       On entry ind specifies the indices of subspace variables.
!       On exit ind is unchanged.
!
!     nenter is an integer variable.
!       On entry nenter is the number of variables entering the
!         free set.
!       On exit nenter is unchanged.
!
!     ileave is an integer variable.
!       On entry indx2(ileave),...,indx2(n) are the variables leaving
!         the free set.
!       On exit ileave is unchanged.
!
!     indx2 is an integer array of dimension n.
!       On entry indx2(1),...,indx2(nenter) are the variables entering
!         the free set, while indx2(ileave),...,indx2(n) are the
!         variables leaving the free set.
!       On exit indx2 is unchanged.
!
!     iupdat is an integer variable.
!       On entry iupdat is the total number of BFGS updates made so far.
!       On exit iupdat is unchanged.
!
!     updatd is a logical variable.
!       On entry 'updatd' is true if the L-BFGS matrix is updatd.
!       On exit 'updatd' is unchanged.
!
!     wn is a double precision array of dimension 2m x 2m.
!       On entry wn is unspecified.
!       On exit the upper triangle of wn stores the LEL^T factorization
!         of the 2*col x 2*col indefinite matrix
!                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!
!     wn1 is a double precision array of dimension 2m x 2m.
!       On entry wn1 stores the lower triangular part of
!                     [Y' ZZ'Y   L_a'+R_z']
!                     [L_a+R_z   S'AA'S   ]
!         in the previous iteration.
!       On exit wn1 stores the corresponding updated matrices.
!       The purpose of wn1 is just to store these inner products
!       so they can be easily updated and inserted into wn.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     ws, wy, sy, and wtyy are double precision arrays;
!     theta is a double precision variable;
!     col is an integer variable;
!     head is an integer variable.
!       On entry they store the information defining the
!                                          limited memory BFGS matrix:
!         ws(n,m) stores S, a set of s-vectors;
!         wy(n,m) stores Y, a set of y-vectors;
!         sy(m,m) stores S'Y;
!         wtyy(m,m) stores the Cholesky factorization
!                                   of (theta*S'S+LD^(-1)L')
!         theta is the scaling factor specifying B_0 = theta I;
!         col is the number of variable metric corrections stored;
!         head is the location of the 1st s- (or y-) vector in S (or Y).
!       On exit they are unchanged.
!
!     info is an integer variable.
!       On entry info is unspecified.
!       On exit info =  0 for normal return;
!                    = -1 when the 1st Cholesky factorization failed;
!                    = -2 when the 2st Cholesky factorization failed.
!
!     Subprograms called:
!
!       Linpack ... dcopy, dpofa, dtrsl.
!
!
!     References:
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
!       limited memory FORTRAN code for solving bound constrained
!       optimization problems'', Tech. Report, NAM-11, EECS Department,
!       Northwestern University, 1994.
!
!       (Postscript files of these papers are available via anonymous
!        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!     Form the lower triangular part of
!               WN1 = [Y' ZZ'Y   L_a'+R_z']
!                     [L_a+R_z   S'AA'S   ]
!        where L_a is the strictly lower triangular part of S'AA'Y
!              R_z is the upper triangular part of S'ZZ'Y.

      IF (updatd) THEN
         IF (iupdat > m) THEN
!                                 shift old part of WN1.
            DO 10 jy = 1, m - 1
               js = m + jy
               CALL dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
               CALL dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
               CALL dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
   10       CONTINUE
         END IF

!          put new rows in blocks (1,1), (2,1) and (2,2).
         pbegin = 1
         pend = nsub
         dbegin = nsub + 1
         dend = n
         iy = col
         is = m + col
         ipntr = head + col - 1
         IF (ipntr > m) ipntr = ipntr - m
         jpntr = head
         DO 20 jy = 1, col
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
!             compute element jy of row 'col' of Y'ZZ'Y
            DO 15 k = pbegin, pend
               k1 = ind(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
   15       CONTINUE
!             compute elements jy of row 'col' of L_a and S'AA'S
            DO 16 k = dbegin, dend
               k1 = ind(k)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
   16       CONTINUE
            wn1(iy,jy) = temp1
            wn1(is,js) = temp2
            wn1(is,jy) = temp3
            jpntr = MOD(jpntr,m) + 1
   20    CONTINUE

!          put new column in block (2,1).
         jy = col
         jpntr = head + col - 1
         IF (jpntr > m) jpntr = jpntr - m
         ipntr = head
         DO 30 i = 1, col
            is = m + i
            temp3 = zero
!             compute element i of column 'col' of R_z
            DO 25 k = pbegin, pend
               k1 = ind(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
   25       CONTINUE
            ipntr = MOD(ipntr,m) + 1
            wn1(is,jy) = temp3
   30    CONTINUE
         upcl = col - 1
      ELSE
         upcl = col
      END IF

!       modify the old parts in blocks (1,1) and (2,2) due to changes
!       in the set of free variables.
      ipntr = head
      DO 45 iy = 1, upcl
         is = m + iy
         jpntr = head
         DO 40 jy = 1, iy
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            temp4 = zero
            DO 35 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
   35       CONTINUE
            DO 36 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
               temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
   36       CONTINUE
            wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
            wn1(is,js) = wn1(is,js) - temp2 + temp4
            jpntr = MOD(jpntr,m) + 1
   40    CONTINUE
         ipntr = MOD(ipntr,m) + 1
   45 CONTINUE

!       modify the old parts in block (2,1).
      ipntr = head
      DO 60 is = m + 1, m + upcl
         jpntr = head
         DO 55 jy = 1, upcl
            temp1 = zero
            temp3 = zero
            DO 50 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
   50       CONTINUE
            DO 51 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
   51       CONTINUE
         IF (is <= jy + m) THEN
               wn1(is,jy) = wn1(is,jy) + temp1 - temp3
            ELSE
               wn1(is,jy) = wn1(is,jy) - temp1 + temp3
            END IF
            jpntr = MOD(jpntr,m) + 1
   55    CONTINUE
         ipntr = MOD(ipntr,m) + 1
   60 CONTINUE

!     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
!                                     [-L_a +R_z        S'AA'S*theta]

      m2 = 2*m
      DO 70 iy = 1, col
         is = col + iy
         is1 = m + iy
         DO 65 jy = 1, iy
            js = col + jy
            js1 = m + jy
            wn(jy,iy) = wn1(iy,jy)/theta
            wn(js,is) = wn1(is1,js1)*theta
   65    CONTINUE
         DO 66 jy = 1, iy - 1
            wn(jy,is) = -wn1(is1,jy)
   66    CONTINUE
         DO 67 jy = iy, col
            wn(jy,is) = wn1(is1,jy)
   67    CONTINUE
         wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
   70 CONTINUE

!     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
!                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

!        first Cholesky factor (1,1) block of wn to get LL'
!                          with L' stored in the upper triangle of wn.
      CALL dpofa(wn,m2,col,info)
      IF (info /= 0) THEN
         info = -1
         RETURN
      END IF
!        then form L^-1(-L_a'+R_z') in the (1,2) block.
      col2 = 2*col
      DO 71 js = col+1 ,col2
         CALL dtrsl(wn,m2,col,wn(1,js),11,info)
   71 CONTINUE

!     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
!        upper triangle of (2,2) block of wn.


      DO 72 is = col+1, col2
         DO 74 js = is, col2
               wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
   74       CONTINUE
   72    CONTINUE

!     Cholesky factorization of (2,2) block of wn.

      CALL dpofa(wn(col+1,col+1),m2,col,info)
      IF (info /= 0) THEN
         info = -2
         RETURN
      END IF

      RETURN

      END SUBROUTINE formk

!======================= The end of formk ==============================

! *****************************************************************************
!> \brief ...
!> \param m ...
!> \param wt ...
!> \param sy ...
!> \param ss ...
!> \param col ...
!> \param theta ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE formt(m, wt, sy, ss, col, theta, info)

    INTEGER                                  :: m
    REAL(KIND=dp)                            :: wt(m, m), sy(m, m), ss(m, m)
    INTEGER                                  :: col
    REAL(KIND=dp)                            :: theta
    INTEGER                                  :: info

    REAL(KIND=dp), PARAMETER                 :: zero = 0.0_dp

    INTEGER                                  :: i, j, k, k1
    REAL(KIND=dp)                            :: ddum

!     ************
!
!     Subroutine formt
!
!       This subroutine forms the upper half of the pos. def. and symm.
!         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
!         of the array wt, and performs the Cholesky factorization of T
!         to produce J*J', with J' stored in the upper triangle of wt.
!
!     Subprograms called:
!
!       Linpack ... dpofa.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!     Form the upper half of  T = theta*SS + L*D^(-1)*L',
!        store T in the upper triangle of the array wt.

      DO 52 j = 1, col
         wt(1,j) = theta*ss(1,j)
   52 CONTINUE
      DO 55 i = 2, col
         DO 54 j = i, col
            k1 = MIN(i,j) - 1
            ddum  = zero
            DO 53 k = 1, k1
               ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
   53       CONTINUE
            wt(i,j) = ddum + theta*ss(i,j)
   54    CONTINUE
   55 CONTINUE

!     Cholesky factorize T to J*J' with
!        J' stored in the upper triangle of wt.

      CALL dpofa(wt,m,col,info)
      IF (info /= 0) THEN
         info = -3
      END IF

      RETURN

      END SUBROUTINE formt

!======================= The end of formt ==============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param nfree ...
!> \param index ...
!> \param nenter ...
!> \param ileave ...
!> \param indx2 ...
!> \param iwhere ...
!> \param wrk ...
!> \param updatd ...
!> \param cnstnd ...
!> \param iprint ...
!> \param iter ...
! *****************************************************************************
      SUBROUTINE freev(n, nfree, index, nenter, ileave, indx2,          &
                       iwhere, wrk, updatd, cnstnd, iprint, iter)

    INTEGER                                  :: n, nfree, INDEX(n), nenter, &
                                                ileave, indx2(n), iwhere(n)
    LOGICAL                                  :: wrk, updatd, cnstnd
    INTEGER                                  :: iprint, iter

    INTEGER                                  :: i, iact, k

!     ************
!
!     Subroutine freev
!
!     This subroutine counts the entering and leaving variables when
!       iter > 0, and finds the index set of free and active variables
!       at the GCP.
!
!     cnstnd is a logical variable indicating whether bounds are present
!
!     index is an integer array of dimension n
!       for i=1,...,nfree, index(i) are the indices of free variables
!       for i=nfree+1,...,n, index(i) are the indices of bound variables
!       On entry after the first iteration, index gives
!         the free variables at the previous iteration.
!       On exit it gives the free variables based on the determination
!         in cauchy using the array iwhere.
!
!     indx2 is an integer array of dimension n
!       On entry indx2 is unspecified.
!       On exit with iter>0, indx2 indicates which variables
!          have changed status since the previous iteration.
!       For i= 1,...,nenter, indx2(i) have changed from bound to free.
!       For i= ileave+1,...,n, indx2(i) have changed from free to bound.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      nenter = 0
      ileave = n + 1
      IF (iter > 0 .AND. cnstnd) THEN
!                           count the entering and leaving variables.
         DO 20 i = 1, nfree
            k = INDEX(i)

!            write(*,*) ' k  = index(i) ', k
!            write(*,*) ' index = ', i

            IF (iwhere(k) > 0) THEN
               ileave = ileave - 1
               indx2(ileave) = k
               IF (iprint >= 100) WRITE (*,*)                         &
     &             'Variable ',k,' leaves the set of free variables'
            END IF
   20    CONTINUE
         DO 22 i = 1 + nfree, n
            k = INDEX(i)
            IF (iwhere(k) <= 0) THEN
               nenter = nenter + 1
               indx2(nenter) = k
               IF (iprint >= 100) WRITE (*,*)                         &
     &             'Variable ',k,' enters the set of free variables'
            END IF
   22    CONTINUE
         IF (iprint >= 99) WRITE (*,*)                                &
             n+1-ileave,' variables leave; ',nenter,' variables enter'
      END IF
      wrk = (ileave < n+1) .OR. (nenter > 0) .OR. updatd

!     Find the index set of free and active variables at the GCP.

      nfree = 0
      iact = n + 1
      DO 24 i = 1, n
         IF (iwhere(i) <= 0) THEN
            nfree = nfree + 1
            INDEX(nfree) = i
         ELSE
            iact = iact - 1
            INDEX(iact) = i
         END IF
   24 CONTINUE
      IF (iprint >= 99) WRITE (*,*)                                   &
            nfree,' variables are free at GCP ',iter + 1

      RETURN

      END SUBROUTINE freev

!======================= The end of freev ==============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param t ...
!> \param iorder ...
!> \param iheap ...
! *****************************************************************************
      SUBROUTINE hpsolb(n, t, iorder, iheap)
    INTEGER                                  :: n
    REAL(KIND=dp)                            :: t(n)
    INTEGER                                  :: iorder(n), iheap

    INTEGER                                  :: i, indxin, indxou, j, k
    REAL(KIND=dp)                            :: ddum, out

!     ************
!
!     Subroutine hpsolb
!
!     This subroutine sorts out the least element of t, and puts the
!       remaining elements of t in a heap.
!
!     n is an integer variable.
!       On entry n is the dimension of the arrays t and iorder.
!       On exit n is unchanged.
!
!     t is a double precision array of dimension n.
!       On entry t stores the elements to be sorted,
!       On exit t(n) stores the least elements of t, and t(1) to t(n-1)
!         stores the remaining elements in the form of a heap.
!
!     iorder is an integer array of dimension n.
!       On entry iorder(i) is the index of t(i).
!       On exit iorder(i) is still the index of t(i), but iorder may be
!         permuted in accordance with t.
!
!     iheap is an integer variable specifying the task.
!       On entry iheap should be set as follows:
!         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
!         iheap .ne. 0 if otherwise.
!       On exit iheap is unchanged.
!
!
!     References:
!       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!     ************

      IF (iheap == 0) THEN

!        Rearrange the elements t(1) to t(n) to form a heap.

         DO 20 k = 2, n
            ddum  = t(k)
            indxin = iorder(k)

!           Add ddum to the heap.
            i = k
   10       CONTINUE
            IF (i.gt.1) THEN
               j = i/2
               IF (ddum < t(j)) THEN
                  t(i) = t(j)
                  iorder(i) = iorder(j)
                  i = j
                  GOTO 10
               END IF
            END IF
            t(i) = ddum
            iorder(i) = indxin
   20    CONTINUE
      END IF

!     Assign to 'out' the value of t(1), the least member of the heap,
!        and rearrange the remaining members to form a heap as
!        elements 1 to n-1 of t.

      IF (n > 1) THEN
         i = 1
         out = t(1)
         indxou = iorder(1)
         ddum  = t(n)
         indxin  = iorder(n)

!        Restore the heap
   30    CONTINUE
         j = i+i
         IF (j <= n-1) THEN
            IF (t(j+1) < t(j)) j = j+1
            IF (t(j) < ddum ) THEN
               t(i) = t(j)
               iorder(i) = iorder(j)
               i = j
               GOTO 30
            END IF
         END IF
         t(i) = ddum
         iorder(i) = indxin

!     Put the least member in t(n).

         t(n) = out
         iorder(n) = indxou
      END IF

      RETURN

      END SUBROUTINE hpsolb

!====================== The end of hpsolb ==============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param x ...
!> \param f ...
!> \param fold ...
!> \param gd ...
!> \param gdold ...
!> \param g ...
!> \param d ...
!> \param r ...
!> \param t ...
!> \param z ...
!> \param stp ...
!> \param dnorm ...
!> \param dtd ...
!> \param xstep ...
!> \param stpmx ...
!> \param iter ...
!> \param ifun ...
!> \param iback ...
!> \param nfgv ...
!> \param info ...
!> \param task ...
!> \param boxed ...
!> \param cnstnd ...
!> \param csave ...
!> \param isave ...
!> \param dsave ...
! *****************************************************************************
      SUBROUTINE lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t,&
                        z, stp, dnorm, dtd, xstep, stpmx, iter, ifun,   &
                        iback, nfgv, info, task, boxed, cnstnd, csave,  &
                        isave, dsave)

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp)                            :: x(n), f, fold, gd, gdold, &
                                                g(n), d(n), r(n), t(n), z(n), &
                                                stp, dnorm, dtd, xstep, stpmx
    INTEGER                                  :: iter, ifun, iback, nfgv, info
    CHARACTER(LEN=60)                        :: task
    LOGICAL                                  :: boxed, cnstnd
    CHARACTER(LEN=60)                        :: csave
    INTEGER                                  :: isave(2)
    REAL(KIND=dp)                            :: dsave(13)

    REAL(KIND=dp), PARAMETER :: big = 1.0E10_dp, ftol = 1.0E-3_dp, &
      gtol = 0.9_dp, one = 1.0_dp, xtol = 0.1_dp, zero = 0.0_dp

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: a1, a2, ddot

!     **********
!
!     Subroutine lnsrlb
!
!     This subroutine calls subroutine dcsrch from the Minpack2 library
!       to perform the line search.  Subroutine dscrch is safeguarded so
!       that all trial points lie within the feasible region.
!
!     Subprograms called:
!
!       Minpack2 Library ... dcsrch.
!
!       Linpack ... dtrsl, ddot.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     **********

      IF (task(1:5) == 'FG_LN') GOTO 556

      dtd = ddot(n,d,1,d,1)
      dnorm = SQRT(dtd)

!     Determine the maximum step length.

      stpmx = big
      IF (cnstnd) THEN
         IF (iter == 0) THEN
            stpmx = one
         ELSE
            DO 43 i = 1, n
               a1 = d(i)
               IF (nbd(i) /= 0) THEN
                  IF (a1 < zero .AND. nbd(i) <= 2) THEN
                     a2 = l(i) - x(i)
                     IF (a2 >= zero) THEN
                        stpmx = zero
                     ELSE IF (a1*stpmx < a2) THEN
                        stpmx = a2/a1
                     END IF
                  ELSE IF (a1 > zero .AND. nbd(i) >= 2) THEN
                     a2 = u(i) - x(i)
                     IF (a2 <= zero) THEN
                        stpmx = zero
                     ELSE IF (a1*stpmx > a2) THEN
                        stpmx = a2/a1
                     END IF
                  END IF
               END IF
   43       CONTINUE
         END IF
      END IF

      IF (iter == 0 .AND. .NOT. boxed) THEN
         stp = MIN(one/dnorm, stpmx)
      ELSE
         stp = one
      END IF

      CALL dcopy(n,x,1,t,1)
      CALL dcopy(n,g,1,r,1)
      fold = f
      ifun = 0
      iback = 0
      csave = 'START'
  556 CONTINUE
      gd = ddot(n,g,1,d,1)
      IF (ifun == 0) THEN
         gdold=gd
         IF (gd >= zero) THEN
!                               the directional derivative >=0.
!                               Line search is impossible.
            WRITE(*,*)' ascent direction in projection gd = ', gd
            info = -4
            RETURN
         END IF
      END IF

      CALL dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)

      xstep = stp*dnorm
      IF (csave(1:4) /= 'CONV' .AND. csave(1:4) /= 'WARN') THEN
         task = 'FG_LNSRCH'
         ifun = ifun + 1
         nfgv = nfgv + 1
         iback = ifun - 1
         IF (stp == one) THEN
            CALL dcopy(n,z,1,x,1)
         ELSE
            DO 41 i = 1, n
               x(i) = stp*d(i) + t(i)
   41       CONTINUE
         END IF
      ELSE
         task = 'NEW_X'
      END IF

      RETURN

      END SUBROUTINE lnsrlb

!======================= The end of lnsrlb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param ws ...
!> \param wy ...
!> \param sy ...
!> \param ss ...
!> \param d ...
!> \param r ...
!> \param itail ...
!> \param iupdat ...
!> \param col ...
!> \param head ...
!> \param theta ...
!> \param rr ...
!> \param dr ...
!> \param stp ...
!> \param dtd ...
! *****************************************************************************
      SUBROUTINE matupd(n, m, ws, wy, sy, ss, d, r, itail,              &
                        iupdat, col, head, theta, rr, dr, stp, dtd)

    INTEGER                                  :: n, m
    REAL(KIND=dp)                            :: ws(n, m), wy(n, m), sy(m, m), &
                                                ss(m, m), d(n), r(n)
    INTEGER                                  :: itail, iupdat, col, head
    REAL(KIND=dp)                            :: theta, rr, dr, stp, dtd

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp

    INTEGER                                  :: j, pointr
    REAL(KIND=dp)                            :: ddot

!     ************
!
!     Subroutine matupd
!
!       This subroutine updates matrices WS and WY, and forms the
!         middle matrix in B.
!
!     Subprograms called:
!
!       Linpack ... dcopy, ddot.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!     Set pointers for matrices WS and WY.

      IF (iupdat <= m) THEN
         col = iupdat
         itail = MOD(head+iupdat-2,m) + 1
      ELSE
         itail = MOD(itail,m) + 1
         head = MOD(head,m) + 1
      END IF

!     Update matrices WS and WY.

      CALL dcopy(n,d,1,ws(1,itail),1)
      CALL dcopy(n,r,1,wy(1,itail),1)

!     Set theta=yy/ys.

      theta = rr/dr

!     Form the middle matrix in B.

!        update the upper triangle of SS,
!                                         and the lower triangle of SY:
      IF (iupdat > m) THEN
!                              move old information
         DO 50 j = 1, col - 1
            CALL dcopy(j,ss(2,j+1),1,ss(1,j),1)
            CALL dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
   50    CONTINUE
      END IF
!        add new information: the last row of SY
!                                             and the last column of SS:
      pointr = head
      DO 51 j = 1, col - 1
         sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
         ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
         pointr = MOD(pointr,m) + 1
   51 CONTINUE
      IF (stp == one) THEN
         ss(col,col) = dtd
      ELSE
         ss(col,col) = stp*stp*dtd
      END IF
      sy(col,col) = dr

      RETURN

      END SUBROUTINE matupd

!======================= The end of matupd =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param l ...
!> \param u ...
!> \param x ...
!> \param iprint ...
!> \param itfile ...
!> \param epsmch ...
! *****************************************************************************
      SUBROUTINE prn1lb(n, m, l, u, x, iprint, itfile, epsmch)

    INTEGER                                  :: n, m
    REAL(KIND=dp)                            :: l(n), u(n), x(n)
    INTEGER                                  :: iprint, itfile
    REAL(KIND=dp)                            :: epsmch

    INTEGER                                  :: i

!     ************
!
!     Subroutine prn1lb
!
!     This subroutine prints the input data, initial point, upper and
!       lower bounds of each variable, machine precision, as well as
!       the headings of the output.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      IF (iprint >= 0) THEN
         WRITE (*,7001) epsmch
         WRITE (*,*) 'N = ',n,'    M = ',m
         IF (iprint >= 1) THEN
            WRITE (itfile,2001) epsmch
            WRITE (itfile,*)'N = ',n,'    M = ',m
            WRITE (itfile,9001)
            IF (iprint > 100) THEN
               WRITE (*,1004) 'L =',(l(i),i = 1,n)
               WRITE (*,1004) 'X0 =',(x(i),i = 1,n)
               WRITE (*,1004) 'U =',(u(i),i = 1,n)
            END IF
         END IF
      END IF

 1004 FORMAT (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 FORMAT ('RUNNING THE L-BFGS-B CODE',/,/,                          &
     & 'it    = iteration number',/,                                    &
     & 'nf    = number of function evaluations',/,                      &
     & 'nseg  = number of segments explored during the Cauchy search',/,&
     & 'nact  = number of active bounds at the generalized Cauchy point'&
     & ,/,                                                              &
     & 'sub   = manner in which the subspace minimization terminated:'  &
     & ,/,'        con = converged, bnd = a bound was reached',/,       &
     & 'itls  = number of iterations performed in the line search',/,   &
     & 'stepl = step length used',/,                                    &
     & 'tstep = norm of the displacement (total step)',/,               &
     & 'projg = norm of the projected gradient',/,                      &
     & 'f     = function value',/,/,                                    &
     & '           * * *',/,/,                                          &
     & 'Machine precision =',1p,d10.3)
 7001 FORMAT ('RUNNING THE L-BFGS-B CODE',/,/,                          &
     & '           * * *',/,/,                                          &
     & 'Machine precision =',1p,d10.3)
 9001 FORMAT (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,'itls', &
              2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

      RETURN

      END SUBROUTINE prn1lb

!======================= The end of prn1lb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
!> \param f ...
!> \param g ...
!> \param iprint ...
!> \param itfile ...
!> \param iter ...
!> \param nfgv ...
!> \param nact ...
!> \param sbgnrm ...
!> \param nseg ...
!> \param word ...
!> \param iword ...
!> \param iback ...
!> \param stp ...
!> \param xstep ...
! *****************************************************************************
      SUBROUTINE prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact,   &
                        sbgnrm, nseg, word, iword, iback, stp, xstep)

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: x(n), f, g(n)
    INTEGER                                  :: iprint, itfile, iter, nfgv, &
                                                nact
    REAL(KIND=dp)                            :: sbgnrm
    INTEGER                                  :: nseg
    CHARACTER(LEN=3)                         :: word
    INTEGER                                  :: iword, iback
    REAL(KIND=dp)                            :: stp, xstep

    INTEGER                                  :: i, imod

!     ************
!
!     Subroutine prn2lb
!
!     This subroutine prints out new information after a successful
!       line search.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!           'word' records the status of subspace solutions.

      IF (iword == 0) THEN
!                            the subspace minimization converged.
         word = 'con'
      ELSE IF (iword == 1) THEN
!                          the subspace minimization stopped at a bound.
         word = 'bnd'
      ELSE IF (iword == 5) THEN
!                             the truncated Newton step has been used.
         word = 'TNT'
      ELSE
         word = '---'
      END IF
      IF (iprint >= 99) THEN
         WRITE (*,*) 'LINE SEARCH',iback,' times; norm of step = ',xstep
         WRITE (*,2001) iter,f,sbgnrm
         IF (iprint > 100) THEN
            WRITE (*,1004) 'X =',(x(i), i = 1, n)
            WRITE (*,1004) 'G =',(g(i), i = 1, n)
         END IF
      ELSE IF (iprint > 0) THEN
         imod = MOD(iter,iprint)
         IF (imod == 0) WRITE (*,2001) iter,f,sbgnrm
      END IF
      IF (iprint >= 1) WRITE (itfile,3001)                            &
                iter,nfgv,nseg,nact,word,iback,stp,xstep,sbgnrm,f

 1004 FORMAT (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 FORMAT                                                            &
     &  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 3001 FORMAT(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3))

      RETURN

      END SUBROUTINE prn2lb

!======================= The end of prn2lb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
!> \param f ...
!> \param task ...
!> \param iprint ...
!> \param info ...
!> \param itfile ...
!> \param iter ...
!> \param nfgv ...
!> \param nintol ...
!> \param nskip ...
!> \param nact ...
!> \param sbgnrm ...
!> \param time ...
!> \param nseg ...
!> \param word ...
!> \param iback ...
!> \param stp ...
!> \param xstep ...
!> \param k ...
!> \param cachyt ...
!> \param sbtime ...
!> \param lnscht ...
! *****************************************************************************
      SUBROUTINE prn3lb(n, x, f, task, iprint, info, itfile,            &
                        iter, nfgv, nintol, nskip, nact, sbgnrm,        &
                        time, nseg, word, iback, stp, xstep, k,         &
                        cachyt, sbtime, lnscht)

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: x(n), f
    CHARACTER(LEN=60)                        :: task
    INTEGER                                  :: iprint, info, itfile, iter, &
                                                nfgv, nintol, nskip, nact
    REAL(KIND=dp)                            :: sbgnrm, time
    INTEGER                                  :: nseg
    CHARACTER(LEN=3)                         :: word
    INTEGER                                  :: iback
    REAL(KIND=dp)                            :: stp, xstep
    INTEGER                                  :: k
    REAL(KIND=dp)                            :: cachyt, sbtime, lnscht

    INTEGER                                  :: i

!     ************
!
!     Subroutine prn3lb
!
!     This subroutine prints out information when either a built-in
!       convergence test is satisfied or when an error message is
!       generated.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      IF (task(1:5) == 'ERROR') GOTO 999

      IF (iprint >= 0) THEN
         WRITE(*,3003)
         WRITE(*,3004)
         WRITE(*,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
         IF (iprint >= 100) THEN
            WRITE (*,1004) 'X =',(x(i),i = 1,n)
         END IF
         IF (iprint >= 1) WRITE (*,*) ' F =',f
      END IF
  999 CONTINUE
      IF (iprint >= 0) THEN
         WRITE (*,3009) task
         IF (info /= 0) THEN
            IF (info == -1) WRITE (*,9011)
            IF (info == -2) WRITE (*,9012)
            IF (info == -3) WRITE (*,9013)
            IF (info == -4) WRITE (*,9014)
            IF (info == -5) WRITE (*,9015)
            IF (info == -6) WRITE (*,*)' Input nbd(',k,') is invalid.'
            IF (info == -7)                                           &
            WRITE (*,*)' l(',k,') > u(',k,').  No feasible solution.'
            IF (info == -8) WRITE (*,9018)
            IF (info == -9) WRITE (*,9019)
         END IF
         IF (iprint >= 1) WRITE (*,3007) cachyt,sbtime,lnscht
         WRITE (*,3008) time
         IF (iprint >= 1) THEN
            IF (info == -4 .OR. info == -9) THEN
               WRITE (itfile,3002)                                      &
                   iter,nfgv,nseg,nact,word,iback,stp,xstep
            END IF
            WRITE (itfile,3009) task
            IF (info /= 0) THEN
               IF (info == -1) WRITE (itfile,9011)
               IF (info == -2) WRITE (itfile,9012)
               IF (info == -3) WRITE (itfile,9013)
               IF (info == -4) WRITE (itfile,9014)
               IF (info == -5) WRITE (itfile,9015)
               IF (info == -8) WRITE (itfile,9018)
               IF (info == -9) WRITE (itfile,9019)
            END IF
            WRITE (itfile,3008) time
         END IF
      END IF

 1004 FORMAT (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 3002 FORMAT(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,'-',10x,'-')
 3003 FORMAT (/,                                                        &
     & '           * * *',/,/,                                          &
     & 'Tit   = total number of iterations',/,                          &
     & 'Tnf   = total number of function evaluations',/,                &
     & 'Tnint = total number of segments explored during',              &
     &           ' Cauchy searches',/,                                  &
     & 'Skip  = number of BFGS updates skipped',/,                      &
     & 'Nact  = number of active bounds at final generalized',          &
     &          ' Cauchy point',/,                                      &
     & 'Projg = norm of the final projected gradient',/,                &
     & 'F     = final function value',/,/,                              &
     & '           * * *')
 3004 FORMAT (/,3x,'N',4x,'Tit',5x,'Tnf',2x,'Tnint',2x,                 &
     &       'Skip',2x,'Nact',5x,'Projg',8x,'F')
 3005 FORMAT (i5,2(1x,i6),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
 3007 FORMAT (/,' Cauchy                time',1p,e10.3,' seconds.',/    &
     &        ' Subspace minimization time',1p,e10.3,' seconds.',/      &
     &        ' Line search           time',1p,e10.3,' seconds.')
 3008 FORMAT (/,' Total User time',1p,e10.3,' seconds.',/)
 3009 FORMAT (/,a60)
 9011 FORMAT (/,                                                        &
     &' Matrix in 1st Cholesky factorization in formk is not Pos. Def.')
 9012 FORMAT (/,                                                        &
     &' Matrix in 2st Cholesky factorization in formk is not Pos. Def.')
 9013 FORMAT (/,                                                        &
     &' Matrix in the Cholesky factorization in formt is not Pos. Def.')
 9014 FORMAT (/,                                                        &
     &' Derivative >= 0, backtracking line search impossible.',/,       &
     &'   Previous x, f and g restored.',/,                             &
     &' Possible causes: 1 error in function or gradient evaluation;',/,&
     &'                  2 rounding errors dominate computation.')
 9015 FORMAT (/,                                                        &
     &' Warning:  more than 10 function and gradient',/,                &
     &'   evaluations in the last line search.  Termination',/,         &
     &'   may possibly be caused by a bad search direction.')
 9018 FORMAT (/,' The triangular system is singular.')
 9019 FORMAT (/,                                                        &
     &' Line search cannot locate an adequate point after 20 function',/&
     &,'  and gradient evaluations.  Previous x, f and g restored.',/,  &
     &' Possible causes: 1 error in function or gradient evaluation;',/,&
     &'                  2 rounding error dominate computation.')

      RETURN

      END SUBROUTINE prn3lb

!======================= The end of prn3lb =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param x ...
!> \param g ...
!> \param sbgnrm ...
! *****************************************************************************
      SUBROUTINE projgr(n, l, u, nbd, x, g, sbgnrm)

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp)                            :: x(n), g(n), sbgnrm

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: gi

!     ************
!
!     Subroutine projgr
!
!     This subroutine computes the infinity norm of the projected
!       gradient.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      sbgnrm = zero
      DO 15 i = 1, n
        gi = g(i)
        IF (nbd(i) /= 0) THEN
           IF (gi < zero) THEN
              IF (nbd(i) >= 2) gi = MAX((x(i)-u(i)),gi)
           ELSE
              IF (nbd(i) <= 2) gi = MIN((x(i)-l(i)),gi)
           END IF
        END IF
        sbgnrm = MAX(sbgnrm,ABS(gi))
   15 CONTINUE

      RETURN

      END SUBROUTINE projgr

!======================= The end of projgr =============================

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param m ...
!> \param nsub ...
!> \param ind ...
!> \param l ...
!> \param u ...
!> \param nbd ...
!> \param x ...
!> \param d ...
!> \param xp ...
!> \param ws ...
!> \param wy ...
!> \param theta ...
!> \param xx ...
!> \param gg ...
!> \param col ...
!> \param head ...
!> \param iword ...
!> \param wv ...
!> \param wn ...
!> \param iprint ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE subsm ( n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy,  &
                         theta, xx, gg,                                 &
                         col, head, iword, wv, wn, iprint, info )
    INTEGER                                  :: n, m, nsub, ind(nsub)
    REAL(KIND=dp)                            :: l(n), u(n)
    INTEGER                                  :: nbd(n)
    REAL(KIND=dp)                            :: x(n), d(n), xp(n), ws(n, m), &
                                                wy(n, m), theta, xx(n), gg(n)
    INTEGER                                  :: col, head, iword
    REAL(KIND=dp)                            :: wv(2*m), wn(2*m, 2*m)
    INTEGER                                  :: iprint, info

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: col2, i, ibd, j, js, jy, k, &
                                                m2, pointr
    REAL(KIND=dp)                            :: alpha, dd_p, dk, temp1, &
                                                temp2, xk

!     ******************************************************************
!
!     This routine contains the major changes in the updated version.
!     The changes are described in the accompanying paper
!
!      Jose Luis Morales, Jorge Nocedal
!      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large
!       Bound Constrained Optimization". Decemmber 27, 2010.
!
!             J.L. Morales  Departamento de Matematicas,
!                           Instituto Tecnologico Autonomo de Mexico
!                           Mexico D.F.
!
!             J, Nocedal    Department of Electrical Engineering and
!                           Computer Science.
!                           Northwestern University. Evanston, IL. USA
!
!                           January 17, 2011
!
!      *****************************************************************
!
!
!     Subroutine subsm
!
!     Given xcp, l, u, r, an index set that specifies
!       the active set at xcp, and an l-BFGS matrix B
!       (in terms of WY, WS, SY, WT, head, col, and theta),
!       this subroutine computes an approximate solution
!       of the subspace problem
!
!       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
!
!             subject to l<=x<=u
!                       x_i=xcp_i for all i in A(xcp)
!
!       along the subspace unconstrained Newton direction
!
!          d = -(Z'BZ)^(-1) r.
!
!       The formula for the Newton direction, given the L-BFGS matrix
!       and the Sherman-Morrison formula, is
!
!          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
!
!       where
!                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!
!     Note that this procedure for computing d differs
!     from that described in [1]. One can show that the matrix K is
!     equal to the matrix M^[-1]N in that paper.
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     nsub is an integer variable.
!       On entry nsub is the number of free variables.
!       On exit nsub is unchanged.
!
!     ind is an integer array of dimension nsub.
!       On entry ind specifies the coordinate indices of free variables.
!       On exit ind is unchanged.
!
!     l is a double precision array of dimension n.
!       On entry l is the lower bound of x.
!       On exit l is unchanged.
!
!     u is a double precision array of dimension n.
!       On entry u is the upper bound of x.
!       On exit u is unchanged.
!
!     nbd is a integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds, and
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     x is a double precision array of dimension n.
!       On entry x specifies the Cauchy point xcp.
!       On exit x(i) is the minimizer of Q over the subspace of
!                                                        free variables.
!
!     d is a double precision array of dimension n.
!       On entry d is the reduced gradient of Q at xcp.
!       On exit d is the Newton direction of Q.
!
!    xp is a double precision array of dimension n.
!       used to safeguard the projected Newton direction
!
!    xx is a double precision array of dimension n
!       On entry it holds the current iterate
!       On output it is unchanged
!    gg is a double precision array of dimension n
!       On entry it holds the gradient at the current iterate
!       On output it is unchanged
!
!     ws and wy are double precision arrays;
!     theta is a double precision variable;
!     col is an integer variable;
!     head is an integer variable.
!       On entry they store the information defining the
!                                          limited memory BFGS matrix:
!         ws(n,m) stores S, a set of s-vectors;
!         wy(n,m) stores Y, a set of y-vectors;
!         theta is the scaling factor specifying B_0 = theta I;
!         col is the number of variable metric corrections stored;
!         head is the location of the 1st s- (or y-) vector in S (or Y).
!       On exit they are unchanged.
!
!     iword is an integer variable.
!       On entry iword is unspecified.
!       On exit iword specifies the status of the subspace solution.
!         iword = 0 if the solution is in the box,
!                 1 if some bound is encountered.
!
!     wv is a double precision working array of dimension 2m.
!
!     wn is a double precision array of dimension 2m x 2m.
!       On entry the upper triangle of wn stores the LEL^T factorization
!         of the indefinite matrix
!
!              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                  [L_a -R_z           theta*S'AA'S ]
!                                                    where E = [-I  0]
!                                                              [ 0  I]
!       On exit wn is unchanged.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     info is an integer variable.
!       On entry info is unspecified.
!       On exit info = 0       for normal return,
!                    = nonzero for abnormal return
!                                  when the matrix K is ill-conditioned.
!
!     Subprograms called:
!
!       Linpack dtrsl.
!
!
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************
!

      IF (nsub <= 0) RETURN
      IF (iprint >= 99) WRITE (*,1001)

!     Compute wv = W'Zd.

      pointr = head
      DO 20 i = 1, col
         temp1 = zero
         temp2 = zero
         DO 10 j = 1, nsub
            k = ind(j)
            temp1 = temp1 + wy(k,pointr)*d(j)
            temp2 = temp2 + ws(k,pointr)*d(j)
   10    CONTINUE
         wv(i) = temp1
         wv(col + i) = theta*temp2
         pointr = MOD(pointr,m) + 1
   20 CONTINUE

!     Compute wv:=K^(-1)wv.

      m2 = 2*m
      col2 = 2*col
      CALL dtrsl(wn,m2,col2,wv,11,info)
      IF (info /= 0) RETURN
      DO 25 i = 1, col
         wv(i) = -wv(i)
   25    CONTINUE
      CALL dtrsl(wn,m2,col2,wv,01,info)
      IF (info /= 0) RETURN

!     Compute d = (1/theta)d + (1/theta**2)Z'W wv.

      pointr = head
      DO 40 jy = 1, col
         js = col + jy
         DO 30 i = 1, nsub
            k = ind(i)
            d(i) = d(i) + wy(k,pointr)*wv(jy)/theta                     &
     &                  + ws(k,pointr)*wv(js)
   30    CONTINUE
         pointr = MOD(pointr,m) + 1
   40 CONTINUE

      CALL dscal( nsub, one/theta, d, 1 )
!
!-----------------------------------------------------------------
!     Let us try the projection, d is the Newton direction

      iword = 0

      CALL dcopy ( n, x, 1, xp, 1 )
!
      DO 50 i=1, nsub
         k  = ind(i)
         dk = d(i)
         xk = x(k)
         IF ( nbd(k) /= 0 ) THEN
!
                                             ! lower bounds only
            IF ( nbd(k).eq.1 ) THEN
               x(k) = MAX( l(k), xk + dk )
               IF ( x(k).eq.l(k) ) iword = 1
            ELSE
!
                                             ! upper and lower bounds
               IF ( nbd(k).eq.2 ) THEN
                  xk   = MAX( l(k), xk + dk )
                  x(k) = MIN( u(k), xk )
                  IF ( x(k).eq.l(k) .OR. x(k).eq.u(k) ) iword = 1
               ELSE
!
                                             ! upper bounds only
                  IF ( nbd(k).eq.3 ) THEN
                     x(k) = MIN( u(k), xk + dk )
                     IF ( x(k).eq.u(k) ) iword = 1
                  END IF
               END IF
            END IF
!
                                             ! free variables
         ELSE
            x(k) = xk + dk
         END IF
   50 CONTINUE
!
      IF ( iword.eq.0 ) THEN
         go to 911
      END IF
!
!     check sign of the directional derivative
!
      dd_p = zero
      DO 55 i=1, n
         dd_p  = dd_p + (x(i) - xx(i))*gg(i)
   55 CONTINUE
      IF ( dd_p .gt.zero ) THEN
         CALL dcopy( n, xp, 1, x, 1 )
         WRITE(*,*) ' Positive dir derivative in projection '
         WRITE(*,*) ' Using the backtracking step '
      ELSE
         go to 911
      END IF
!
!-----------------------------------------------------------------
!
      alpha = one
      temp1 = alpha
      ibd   = 0
      DO 60 i = 1, nsub
         k = ind(i)
         dk = d(i)
         IF (nbd(k) /= 0) THEN
            IF (dk < zero .AND. nbd(k) <= 2) THEN
               temp2 = l(k) - x(k)
               IF (temp2 >= zero) THEN
                  temp1 = zero
               ELSE IF (dk*alpha < temp2) THEN
                  temp1 = temp2/dk
               END IF
            ELSE IF (dk > zero .AND. nbd(k) >= 2) THEN
               temp2 = u(k) - x(k)
               IF (temp2 <= zero) THEN
                  temp1 = zero
               ELSE IF (dk*alpha > temp2) THEN
                  temp1 = temp2/dk
               END IF
            END IF
            IF (temp1 < alpha) THEN
               alpha = temp1
               ibd = i
            END IF
         END IF
   60 CONTINUE

      IF (alpha < one) THEN
         dk = d(ibd)
         k = ind(ibd)
         IF (dk > zero) THEN
            x(k) = u(k)
            d(ibd) = zero
         ELSE IF (dk < zero) THEN
            x(k) = l(k)
            d(ibd) = zero
         END IF
      END IF
      DO 70 i = 1, nsub
         k    = ind(i)
         x(k) = x(k) + alpha*d(i)
   70 CONTINUE
!ccccc
  911 CONTINUE

      IF (iprint >= 99) WRITE (*,1004)

 1001 FORMAT (/,'----------------SUBSM entered-----------------',/)
 1004 FORMAT (/,'----------------exit SUBSM --------------------',/)

      RETURN

      END SUBROUTINE subsm
!====================== The end of subsm ===============================

! *****************************************************************************
!> \brief ...
!> \param f ...
!> \param g ...
!> \param stp ...
!> \param ftol ...
!> \param gtol ...
!> \param xtol ...
!> \param stpmin ...
!> \param stpmax ...
!> \param task ...
!> \param isave ...
!> \param dsave ...
! *****************************************************************************
      SUBROUTINE dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,           &
                        task,isave,dsave)
    REAL(KIND=dp)                            :: f, g, stp, ftol, gtol, xtol, &
                                                stpmin, stpmax
    CHARACTER(LEN=*)                         :: task
    INTEGER                                  :: isave(2)
    REAL(KIND=dp)                            :: dsave(13)

    REAL(KIND=dp), PARAMETER                 :: p5 = 0.5_dp, p66 = 0.66_dp, &
                                                xtrapl = 1.1_dp, &
                                                xtrapu = 4.0_dp, zero = 0.0_dp

    INTEGER                                  :: stage
    LOGICAL                                  :: brackt
    REAL(KIND=dp)                            :: finit, fm, ftest, fx, fxm, &
                                                fy, fym, ginit, gm, gtest, &
                                                gx, gxm, gy, gym, stmax, &
                                                stmin, stx, sty, width, width1

!     **********
!
!     Subroutine dcsrch
!
!     This subroutine finds a step that satisfies a sufficient
!     decrease condition and a curvature condition.
!
!     Each call of the subroutine updates an interval with
!     endpoints stx and sty. The interval is initially chosen
!     so that it contains a minimizer of the modified function
!
!           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
!
!     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
!     interval is chosen so that it contains a minimizer of f.
!
!     The algorithm is designed to find a step that satisfies
!     the sufficient decrease condition
!
!           f(stp) <= f(0) + ftol*stp*f'(0),
!
!     and the curvature condition
!
!           abs(f'(stp)) <= gtol*abs(f'(0)).
!
!     If ftol is less than gtol and if, for example, the function
!     is bounded below, then there is always a step which satisfies
!     both conditions.
!
!     If no step can be found that satisfies both conditions, then
!     the algorithm stops with a warning. In this case stp only
!     satisfies the sufficient decrease condition.
!
!     A typical invocation of dcsrch has the following outline:
!
!     task = 'START'
!  10 continue
!        call dcsrch( ... )
!        if (task .eq. 'FG') then
!           Evaluate the function and the gradient at stp
!           goto 10
!           end if
!
!     NOTE: The user must no alter work arrays between calls.
!
!     The subroutine statement is
!
!        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
!                          task,isave,dsave)
!     where
!
!       f is a double precision variable.
!         On initial entry f is the value of the function at 0.
!            On subsequent entries f is the value of the
!            function at stp.
!         On exit f is the value of the function at stp.
!
!       g is a double precision variable.
!         On initial entry g is the derivative of the function at 0.
!            On subsequent entries g is the derivative of the
!            function at stp.
!         On exit g is the derivative of the function at stp.
!
!       stp is a double precision variable.
!         On entry stp is the current estimate of a satisfactory
!            step. On initial entry, a positive initial estimate
!            must be provided.
!         On exit stp is the current estimate of a satisfactory step
!            if task = 'FG'. If task = 'CONV' then stp satisfies
!            the sufficient decrease and curvature condition.
!
!       ftol is a double precision variable.
!         On entry ftol specifies a nonnegative tolerance for the
!            sufficient decrease condition.
!         On exit ftol is unchanged.
!
!       gtol is a double precision variable.
!         On entry gtol specifies a nonnegative tolerance for the
!            curvature condition.
!         On exit gtol is unchanged.
!
!       xtol is a double precision variable.
!         On entry xtol specifies a nonnegative relative tolerance
!            for an acceptable step. The subroutine exits with a
!            warning if the relative difference between sty and stx
!            is less than xtol.
!         On exit xtol is unchanged.
!
!       stpmin is a double precision variable.
!         On entry stpmin is a nonnegative lower bound for the step.
!         On exit stpmin is unchanged.
!
!       stpmax is a double precision variable.
!         On entry stpmax is a nonnegative upper bound for the step.
!         On exit stpmax is unchanged.
!
!       task is a character variable of length at least 60.
!         On initial entry task must be set to 'START'.
!         On exit task indicates the required action:
!
!            If task(1:2) = 'FG' then evaluate the function and
!            derivative at stp and call dcsrch again.
!
!            If task(1:4) = 'CONV' then the search is successful.
!
!            If task(1:4) = 'WARN' then the subroutine is not able
!            to satisfy the convergence conditions. The exit value of
!            stp contains the best point found during the search.
!
!            If task(1:5) = 'ERROR' then there is an error in the
!            input arguments.
!
!         On exit with convergence, a warning or an error, the
!            variable task contains additional information.
!
!       isave is an integer work array of dimension 2.
!
!       dsave is a double precision work array of dimension 13.
!
!     Subprograms called
!
!       MINPACK-2 ... dcstep
!
!     MINPACK-1 Project. June 1983.
!     Argonne National Laboratory.
!     Jorge J. More' and David J. Thuente.
!
!     MINPACK-2 Project. October 1993.
!     Argonne National Laboratory and University of Minnesota.
!     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
!
!     **********
!     Initialization block.

      IF (task(1:5) == 'START') THEN

!        Check the input arguments for errors.

         IF (stp < stpmin) task = 'ERROR: STP < STPMIN'
         IF (stp > stpmax) task = 'ERROR: STP > STPMAX'
         IF (g >= zero) task = 'ERROR: INITIAL G >= ZERO'
         IF (ftol < zero) task = 'ERROR: FTOL < ZERO'
         IF (gtol < zero) task = 'ERROR: GTOL < ZERO'
         IF (xtol < zero) task = 'ERROR: XTOL < ZERO'
         IF (stpmin < zero) task = 'ERROR: STPMIN < ZERO'
         IF (stpmax < stpmin) task = 'ERROR: STPMAX < STPMIN'

!        Exit if there are errors on input.

         IF (task(1:5) == 'ERROR') RETURN

!        Initialize local variables.

         brackt = .FALSE.
         stage = 1
         finit = f
         ginit = g
         gtest = ftol*ginit
         width = stpmax - stpmin
         width1 = width/p5

!        The variables stx, fx, gx contain the values of the step,
!        function, and derivative at the best step.
!        The variables sty, fy, gy contain the value of the step,
!        function, and derivative at sty.
!        The variables stp, f, g contain the values of the step,
!        function, and derivative at stp.

         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = stp + xtrapu*stp
         task = 'FG'

         GOTO 1000

      ELSE

!        Restore local variables.

         IF (isave(1) == 1) THEN
            brackt = .TRUE.
         ELSE
            brackt = .FALSE.
         END IF
         stage = isave(2)
         ginit = dsave(1)
         gtest = dsave(2)
         gx = dsave(3)
         gy = dsave(4)
         finit = dsave(5)
         fx = dsave(6)
         fy = dsave(7)
         stx = dsave(8)
         sty = dsave(9)
         stmin = dsave(10)
         stmax = dsave(11)
         width = dsave(12)
         width1 = dsave(13)

      END IF

!     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
!     algorithm enters the second stage.

      ftest = finit + stp*gtest
      IF (stage == 1 .AND. f <= ftest .AND. g >= zero)            &
         stage = 2

!     Test for warnings.

      IF (brackt .AND. (stp <= stmin .OR. stp >= stmax))            &
         task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      IF (brackt .AND. stmax - stmin <= xtol*stmax)                   &
         task = 'WARNING: XTOL TEST SATISFIED'
      IF (stp == stpmax .AND. f <= ftest .AND. g <= gtest)        &
         task = 'WARNING: STP = STPMAX'
      IF (stp == stpmin .AND. (f > ftest .OR. g >= gtest))       &
         task = 'WARNING: STP = STPMIN'

!     Test for convergence.

      IF (f <= ftest .AND. ABS(g) <= gtol*(-ginit))                 &
         task = 'CONVERGENCE'

!     Test for termination.

      IF (task(1:4) == 'WARN' .OR. task(1:4) == 'CONV') GOTO 1000

!     A modified function is used to predict the step during the
!     first stage if a lower function value has been obtained but
!     the decrease is not sufficient.

      IF (stage == 1 .AND. f <= fx .AND. f > ftest) THEN

!        Define the modified function and derivative values.

         fm = f - stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

!        Call dcstep to update stx, sty, and to compute the new step.

         CALL dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,                 &
                     brackt,stmin,stmax)

!        Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      ELSE

!       Call dcstep to update stx, sty, and to compute the new step.

        CALL dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,                        &
                    brackt,stmin,stmax)

      END IF

!     Decide if a bisection step is needed.

      IF (brackt) THEN
         IF (ABS(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx)
         width1 = width
         width = ABS(sty-stx)
      END IF

!     Set the minimum and maximum steps allowed for stp.

      IF (brackt) THEN
         stmin = MIN(stx,sty)
         stmax = MAX(stx,sty)
      ELSE
         stmin = stp + xtrapl*(stp - stx)
         stmax = stp + xtrapu*(stp - stx)
      END IF

!     Force the step to be within the bounds stpmax and stpmin.

      stp = MAX(stp,stpmin)
      stp = MIN(stp,stpmax)

!     If further progress is not possible, let stp be the best
!     point obtained during the search.

      IF (brackt .AND. (stp <= stmin .OR. stp >= stmax)             &
     &   .OR. (brackt .AND. stmax-stmin <= xtol*stmax)) stp = stx

!     Obtain another function and derivative.

      task = 'FG'

 1000 CONTINUE

!     Save local variables.

      IF (brackt) THEN
         isave(1) = 1
      ELSE
         isave(1) = 0
      END IF
      isave(2) = stage
      dsave(1) =  ginit
      dsave(2) =  gtest
      dsave(3) =  gx
      dsave(4) =  gy
      dsave(5) =  finit
      dsave(6) =  fx
      dsave(7) =  fy
      dsave(8) =  stx
      dsave(9) =  sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1

      RETURN
      END SUBROUTINE dcsrch

!====================== The end of dcsrch ==============================

! *****************************************************************************
!> \brief ...
!> \param stx ...
!> \param fx ...
!> \param dx ...
!> \param sty ...
!> \param fy ...
!> \param dy ...
!> \param stp ...
!> \param fp ...
!> \param dp_loc ...
!> \param brackt ...
!> \param stpmin ...
!> \param stpmax ...
! *****************************************************************************
      SUBROUTINE dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp_loc,brackt,           &
                        stpmin,stpmax)
    REAL(KIND=dp)                            :: stx, fx, dx, sty, fy, dy, &
                                                stp, fp, dp_loc
    LOGICAL                                  :: brackt
    REAL(KIND=dp)                            :: stpmin, stpmax

    REAL(KIND=dp), PARAMETER                 :: p66 = 0.66_dp, &
                                                three = 3.0_dp, two = 2.0_dp, &
                                                zero = 0.0_dp

    REAL(KIND=dp)                            :: gamma, p, q, r, s, sgnd, &
                                                stpc, stpf, stpq, theta

!     **********
!
!     Subroutine dcstep
!
!     This subroutine computes a safeguarded step for a search
!     procedure and updates an interval that contains a step that
!     satisfies a sufficient decrease and a curvature condition.
!
!     The parameter stx contains the step with the least function
!     value. If brackt is set to .true. then a minimizer has
!     been bracketed in an interval with endpoints stx and sty.
!     The parameter stp contains the current step.
!     The subroutine assumes that if brackt is set to .true. then
!
!           min(stx,sty) < stp < max(stx,sty),
!
!     and that the derivative at stx is negative in the direction
!     of the step.
!
!     The subroutine statement is
!
!       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp_loc,brackt,
!                         stpmin,stpmax)
!
!     where
!
!       stx is a double precision variable.
!         On entry stx is the best step obtained so far and is an
!            endpoint of the interval that contains the minimizer.
!         On exit stx is the updated best step.
!
!       fx is a double precision variable.
!         On entry fx is the function at stx.
!         On exit fx is the function at stx.
!
!       dx is a double precision variable.
!         On entry dx is the derivative of the function at
!            stx. The derivative must be negative in the direction of
!            the step, that is, dx and stp - stx must have opposite
!            signs.
!         On exit dx is the derivative of the function at stx.
!
!       sty is a double precision variable.
!         On entry sty is the second endpoint of the interval that
!            contains the minimizer.
!         On exit sty is the updated endpoint of the interval that
!            contains the minimizer.
!
!       fy is a double precision variable.
!         On entry fy is the function at sty.
!         On exit fy is the function at sty.
!
!       dy is a double precision variable.
!         On entry dy is the derivative of the function at sty.
!         On exit dy is the derivative of the function at the exit sty.
!
!       stp is a double precision variable.
!         On entry stp is the current step. If brackt is set to .true.
!            then on input stp must be between stx and sty.
!         On exit stp is a new trial step.
!
!       fp is a double precision variable.
!         On entry fp is the function at stp
!         On exit fp is unchanged.
!
!       dp_loc is a double precision variable.
!         On entry dp_loc is the the derivative of the function at stp.
!         On exit dp_loc is unchanged.
!
!       brackt is an logical variable.
!         On entry brackt specifies if a minimizer has been bracketed.
!            Initially brackt must be set to .false.
!         On exit brackt specifies if a minimizer has been bracketed.
!            When a minimizer is bracketed brackt is set to .true.
!
!       stpmin is a double precision variable.
!         On entry stpmin is a lower bound for the step.
!         On exit stpmin is unchanged.
!
!       stpmax is a double precision variable.
!         On entry stpmax is an upper bound for the step.
!         On exit stpmax is unchanged.
!
!     MINPACK-1 Project. June 1983
!     Argonne National Laboratory.
!     Jorge J. More' and David J. Thuente.
!
!     MINPACK-2 Project. October 1993.
!     Argonne National Laboratory and University of Minnesota.
!     Brett M. Averick and Jorge J. More'.
!
!     **********

      sgnd = dp_loc*(dx/ABS(dx))

!     First case: A higher function value. The minimum is bracketed.
!     If the cubic step is closer to stx than the quadratic step, the
!     cubic step is taken, otherwise the average of the cubic and
!     quadratic steps is taken.

      IF (fp > fx) THEN
         theta = three*(fx - fp)/(stp - stx) + dx + dp_loc
         s = MAX(ABS(theta),ABS(dx),ABS(dp_loc))
         gamma = s*SQRT((theta/s)**2 - (dx/s)*(dp_loc/s))
         IF (stp < stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp_loc
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)*          &
     &                                                       (stp - stx)
         IF (ABS(stpc-stx) < ABS(stpq-stx)) THEN
            stpf = stpc
         ELSE
            stpf = stpc + (stpq - stpc)/two
         END IF
         brackt = .TRUE.

!     Second case: A lower function value and derivatives of opposite
!     sign. The minimum is bracketed. If the cubic step is farther from
!     stp than the secant step, the cubic step is taken, otherwise the
!     secant step is taken.

      ELSE IF (sgnd < zero) THEN
         theta = three*(fx - fp)/(stp - stx) + dx + dp_loc
         s = MAX(ABS(theta),ABS(dx),ABS(dp_loc))
         gamma = s*SQRT((theta/s)**2 - (dx/s)*(dp_loc/s))
         IF (stp > stx) gamma = -gamma
         p = (gamma - dp_loc) + theta
         q = ((gamma - dp_loc) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp_loc/(dp_loc - dx))*(stx - stp)
         IF (ABS(stpc-stp) > ABS(stpq-stp)) THEN
            stpf = stpc
         ELSE
            stpf = stpq
         END IF
         brackt = .TRUE.

!     Third case: A lower function value, derivatives of the same sign,
!     and the magnitude of the derivative decreases.

      ELSE IF (ABS(dp_loc) < ABS(dx)) THEN

!        The cubic step is computed only if the cubic tends to infinity
!        in the direction of the step or if the minimum of the cubic
!        is beyond stp. Otherwise the cubic step is defined to be the
!        secant step.

         theta = three*(fx - fp)/(stp - stx) + dx + dp_loc
         s = MAX(ABS(theta),ABS(dx),ABS(dp_loc))

!        The case gamma = 0 only arises if the cubic does not tend
!        to infinity in the direction of the step.

         gamma = s*SQRT(MAX(zero,(theta/s)**2-(dx/s)*(dp_loc/s)))
         IF (stp > stx) gamma = -gamma
         p = (gamma - dp_loc) + theta
         q = (gamma + (dx - dp_loc)) + gamma
         r = p/q
         IF (r < zero .AND. gamma /= zero) THEN
            stpc = stp + r*(stx - stp)
         ELSE IF (stp > stx) THEN
            stpc = stpmax
         ELSE
            stpc = stpmin
         END IF
         stpq = stp + (dp_loc/(dp_loc - dx))*(stx - stp)

         IF (brackt) THEN

!           A minimizer has been bracketed. If the cubic step is
!           closer to stp than the secant step, the cubic step is
!           taken, otherwise the secant step is taken.

            IF (ABS(stpc-stp) < ABS(stpq-stp)) THEN
               stpf = stpc
            ELSE
               stpf = stpq
            END IF
            IF (stp > stx) THEN
               stpf = MIN(stp+p66*(sty-stp),stpf)
            ELSE
               stpf = MAX(stp+p66*(sty-stp),stpf)
            END IF
         ELSE

!           A minimizer has not been bracketed. If the cubic step is
!           farther from stp than the secant step, the cubic step is
!           taken, otherwise the secant step is taken.

            IF (ABS(stpc-stp) > ABS(stpq-stp)) THEN
               stpf = stpc
            ELSE
               stpf = stpq
            END IF
            stpf = MIN(stpmax,stpf)
            stpf = MAX(stpmin,stpf)
         END IF

!     Fourth case: A lower function value, derivatives of the same sign,
!     and the magnitude of the derivative does not decrease. If the
!     minimum is not bracketed, the step is either stpmin or stpmax,
!     otherwise the cubic step is taken.

      ELSE
         IF (brackt) THEN
            theta = three*(fp - fy)/(sty - stp) + dy + dp_loc
            s = MAX(ABS(theta),ABS(dy),ABS(dp_loc))
            gamma = s*SQRT((theta/s)**2 - (dy/s)*(dp_loc/s))
            IF (stp > sty) gamma = -gamma
            p = (gamma - dp_loc) + theta
            q = ((gamma - dp_loc) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         ELSE IF (stp > stx) THEN
            stpf = stpmax
         ELSE
            stpf = stpmin
         END IF
      END IF

!     Update the interval which contains a minimizer.

      IF (fp > fx) THEN
         sty = stp
         fy = fp
         dy = dp_loc
      ELSE
         IF (sgnd < zero) THEN
            sty = stx
            fy = fx
            dy = dx
         END IF
         stx = stp
         fx = fp
         dx = dp_loc
      END IF

!     Compute the new step.

      stp = stpf

      RETURN
      END SUBROUTINE dcstep

!MK LINPACK

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param lda ...
!> \param n ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE dpofa(a,lda,n,info)
    INTEGER                                  :: lda
    REAL(KIND=dp)                            :: a(lda,*)
    INTEGER                                  :: n, info

    INTEGER                                  :: j, jm1, k
    REAL(KIND=dp)                            :: ddot, s, t

!
!     dpofa factors a double precision symmetric positive definite
!     matrix.
!
!     dpofa is usually called by dpoco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the symmetric matrix to be factored.  only the
!                diagonal and upper triangle are used.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix  r  so that  a = trans(r)*r
!                where  trans(r)  is the transpose.
!                the strict lower triangle is unaltered.
!                if  info .ne. 0 , the factorization is not complete.
!
!        info    integer
!                = 0  for normal return.
!                = k  signals an error condition.  the leading minor
!                     of order  k  is not positive definite.
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas ddot
!     fortran sqrt
!
!     internal variables
!
!     begin block with ...exits to 40
!
!

         DO 30 j = 1, n
            info = j
            s = 0.0_dp
            jm1 = j - 1
            IF (jm1 < 1) go to 20
            DO 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       CONTINUE
   20       CONTINUE
            s = a(j,j) - s
!     ......exit
            IF (s <= 0.0_dp) go to 40
            a(j,j) = SQRT(s)
   30    CONTINUE
         info = 0
   40 CONTINUE
      RETURN
      END SUBROUTINE dpofa

!====================== The end of dpofa ===============================

! *****************************************************************************
!> \brief ...
!> \param t ...
!> \param ldt ...
!> \param n ...
!> \param b ...
!> \param job ...
!> \param info ...
! *****************************************************************************
      SUBROUTINE dtrsl(t,ldt,n,b,job,info)
    INTEGER                                  :: ldt
    REAL(KIND=dp)                            :: t(ldt,*)
    INTEGER                                  :: n
    REAL(KIND=dp)                            :: b(*)
    INTEGER                                  :: job, info

    INTEGER                                  :: CASE, j, jj
    REAL(KIND=dp)                            :: ddot, temp

!
!
!     dtrsl solves systems of the form
!
!                   t * x = b
!     or
!                   trans(t) * x = b
!
!     where t is a triangular matrix of order n. here trans(t)
!     denotes the transpose of the matrix t.
!
!     on entry
!
!         t         double precision(ldt,n)
!                   t contains the matrix of the system. the zero
!                   elements of the matrix are not referenced, and
!                   the corresponding elements of the array can be
!                   used to store other information.
!
!         ldt       integer
!                   ldt is the leading dimension of the array t.
!
!         n         integer
!                   n is the order of the system.
!
!         b         double precision(n).
!                   b contains the right hand side of the system.
!
!         job       integer
!                   job specifies what kind of system is to be solved.
!                   if job is
!
!                        00   solve t*x=b, t lower triangular,
!                        01   solve t*x=b, t upper triangular,
!                        10   solve trans(t)*x=b, t lower triangular,
!                        11   solve trans(t)*x=b, t upper triangular.
!
!     on return
!
!         b         b contains the solution, if info .eq. 0.
!                   otherwise b is unaltered.
!
!         info      integer
!                   info contains zero if the system is nonsingular.
!                   otherwise info contains the index of
!                   the first zero diagonal element of t.
!
!     linpack. this version dated 08/14/78 .
!     g. w. stewart, university of maryland, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!     fortran mod
!
!     internal variables
!
!
!     begin block permitting ...exits to 150
!
!        check for zero diagonal elements.
!

         DO 10 info = 1, n
!     ......exit
            IF (t(info,info) == 0.0_dp) go to 150
   10    CONTINUE
         info = 0
!
!        determine the task and go to it.
!
         CASE = 1
         IF (MOD(job,10) /= 0) CASE = 2
         IF (MOD(job,100)/10 /= 0) CASE = CASE + 2

         SELECT CASE(CASE)
         CASE(1)
!
!        solve t*x=b for t lower triangular
!
            b(1) = b(1)/t(1,1)
            IF (n < 2) go to 40
            DO 30 j = 2, n
               temp = -b(j-1)
               CALL daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       CONTINUE
   40       CONTINUE
         CASE(2)
!
!        solve t*x=b for t upper triangular.
!
            b(n) = b(n)/t(n,n)
            IF (n < 2) go to 70
            DO 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               CALL daxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       CONTINUE
   70       CONTINUE
         CASE(3)
!
!        solve trans(t)*x=b for t lower triangular.
!
            b(n) = b(n)/t(n,n)
            IF (n < 2) go to 100
            DO 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/t(j,j)
   90       CONTINUE
  100       CONTINUE
         CASE(4)
!
!        solve trans(t)*x=b for t upper triangular.
!
  110    CONTINUE
            b(1) = b(1)/t(1,1)
            IF (n < 2) go to 130
            DO 120 j = 2, n
               b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/t(j,j)
  120       CONTINUE
  130       CONTINUE
         CASE DEFAULT
            STOP "dtrsl: unexpected case"
         END SELECT

  150 CONTINUE
      RETURN
      END SUBROUTINE dtrsl

!MK Timer

! *****************************************************************************
!> \brief ...
!> \param ttime ...
! *****************************************************************************
      SUBROUTINE timer(ttime)
    REAL(KIND=dp)                            :: ttime

!
!     REAL temp
!
!     This routine computes cpu time in double precision; it makes use o
!     the intrinsic f90 cpu_time therefore a conversion type is
!     needed.
!
!           J.L Morales  Departamento de Matematicas,
!                        Instituto Tecnologico Autonomo de Mexico
!                        Mexico D.F.
!
!           J.L Nocedal  Department of Electrical Engineering and
!                        Computer Science.
!                        Northwestern University. Evanston, IL. USA
!
!                        January 21, 2011
!
!MK      temp = sngl(ttime)
!MK      CALL cpu_time(temp)
!MK      ttime = DBLE(temp)

         ttime = m_walltime()

      END SUBROUTINE timer

END MODULE cp_lbfgs
