!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE bfgs_optimizer

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_param_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: diamat
  USE message_passing,                 ONLY: mp_bcast
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_parser,                       ONLY: close_file,&
                                             finish_parser,&
                                             open_file,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: geoopt_bfgs

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='bfgs_optimizer'


CONTAINS

!!*** **********************************************************************

  SUBROUTINE geoopt_bfgs(qs_env,globenv,gopt_param)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(gopt_param_type), POINTER           :: gopt_param

    CHARACTER(LEN=*), PARAMETER              :: routine_name = 'geoopt_bfgs'
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    CHARACTER(LEN=80)                        :: title
    INTEGER :: atom, dataunit, handle, hesunit, i, iatom, ikind, indf, istat, &
      it, j, maxiter, natom, ndf, nfree, nkind, output_unit
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: conv, hesrest, ionode
    REAL(KIND=dp)                            :: ediff, emin, eold, etot, &
                                                pred, rad, rat, step
    REAL(KIND=dp), ALLOCATABLE :: dg(:), dr(:), dx(:), eigval(:), &
      eigvec(:,:), g(:), gold(:), hes(:), hes_mat(:,:), work(:), x(:), xold(:)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    ionode      = globenv%ionode
    output_unit = globenv%scr
    maxiter     = gopt_param%max_iter
    hesrest     = gopt_param%hesrest
    conv        = .FALSE.
    rad         = 0.10_dp

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE(gopt_param%to_move(natom))

    gopt_param%to_move = .TRUE.

!   Atom fixing constraints for optimizer

    IF(gopt_param%atmfix)THEN
       CALL get_gopt_atmfix(globenv,natom,nfree,gopt_param%to_move)
    ELSE
      nfree = natom
    END IF

    ndf  = 3 * nfree

    ALLOCATE(x(ndf) )
    ALLOCATE(hes(ndf*(ndf+1)/2))
    ALLOCATE(xold(ndf))
    ALLOCATE(g(ndf))
    ALLOCATE(gold(ndf))
    ALLOCATE(dx(ndf))
    ALLOCATE(hes_mat(ndf,ndf))
    ALLOCATE(dg(ndf))
    ALLOCATE(eigval(ndf))
    ALLOCATE(eigvec(ndf,ndf))
    ALLOCATE(work(ndf))
    ALLOCATE(dr(ndf))

!   *** Store the coordinates ***

    indf = 0

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list)
      natom = SIZE(atom_list)
      DO iatom=1,natom
         atom = atom_list(iatom)
         IF( gopt_param%to_move(particle_set(atom)%iatom) ) THEN
            DO j = 1, 3
              indf = indf + 1
              x(indf) = particle_set(atom)%r(j)
              xold(indf) = x(indf)
            END DO
         END IF
      END DO
    END DO

!  Initialize hessian (hes = unitary matrix ) 
   
   DO indf = 1, ndf
     DO j = 1, indf
        i = indf * (indf - 1)/2 + j
        hes(i) = zero
     END DO
     i = indf * (indf-1)/2 + indf
     hes(i) = one
   END DO 


! Geometry optimization starts now

  IF(ionode)THEN
    WRITE(UNIT=output_unit,FMT="(/,T2,80('*'))")
    WRITE(UNIT=output_unit,FMT="(T2,2('*'),23X,A,23X,2('*'))")&
    "STARTING GEOMETRY OPTIMIZATION"
    WRITE(UNIT=output_unit,FMT="(T2,80('*'))")

    WRITE(UNIT=output_unit,FMT="(/,T2,A,I10)")"No. of free atoms         = ",nfree
    WRITE(UNIT=output_unit,FMT="(T2,A,I10)")  "No. of degrees of freedom = ",ndf

    CALL m_flush(output_unit)
  END IF

  CALL write_qs_particle_coordinates(particle_set,cell,globenv) 
  CALL write_particle_distances(particle_set,cell,globenv)

  IF(ionode)THEN
    CALL open_file(file_name="HESSIAN",&
                   file_status="UNKNOWN",&
                   file_form="UNFORMATTED",&
                   file_action="READWRITE",&
                   unit_number=hesunit)

    CALL open_file(file_name="bfgs.data",&
                   file_status="UNKNOWN",&
                   file_form="FORMATTED",&
                   file_action="WRITE",&
                   unit_number=dataunit)
  END IF

! Calculate Energy & Gradients

  CALL qs_forces(qs_env,globenv)

! Copy forces to g, the gradients

  indf = 0
  DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list)
        natom = SIZE(atom_list)
        DO iatom=1,natom
            atom = atom_list(iatom)
            IF( gopt_param%to_move(particle_set(atom)%iatom) ) THEN
              DO j = 1, 3
                indf = indf + 1
                g(indf) = - particle_set(atom)%f(j)
              END DO
            END IF
        END DO
  END DO

  CALL get_qs_env(qs_env=qs_env,energy=energy)

  etot = energy % total
  emin = energy % total

  it = 0
  
  DO 

     it = it + 1

     IF(it > maxiter )THEN
       IF(ionode)THEN
         WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
         "*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***"
         WRITE(UNIT=output_unit,FMT="(T2,A)")&
         "***        EXITING GEOMETRY OPTIMIZATION         ***" 
         CALL m_flush(output_unit)
       END IF
       EXIT
     END IF

     IF (ionode)THEN
        WRITE(UNIT=output_unit,FMT="(/,T2,26('-'))")
        WRITE(UNIT=output_unit,FMT="(T2,A,I6)")&
       "OPTIMIZATION STEP : ", it
        WRITE(UNIT=output_unit,FMT="(T2,26('-'))")
       CALL m_flush(output_unit)
     END IF


!    ** Hessian update/restarting **

     IF( (it == 1) .AND. hesrest )THEN

       IF(ionode)THEN
         DO i=1,SIZE(hes)
           READ(UNIT=hesunit,IOSTAT=istat)hes(i)
         END DO
         REWIND(hesunit)
         IF(istat /= 0)THEN
           WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
           "**  Error while reading HESSIEN **"
           conv = .TRUE.
         ELSE
           WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
           "** Initial Hessian has been restarted successfully **"
         END IF
  
       END IF

       CALL mp_bcast(hes,globenv%source,globenv%group)

     ELSE IF( it > 1 ) THEN

        DO indf = 1, ndf
          dx(indf) = x(indf) - xold(indf)
          dg(indf) = g(indf) - gold(indf)     
        END DO

        CALL bfgs(ndf,dx,dg,hes,work)

        IF(ionode)THEN
           DO i = 1, SIZE(hes)
             WRITE(UNIT=hesunit)hes(i)
            CALL m_flush(hesunit)
           END DO
           REWIND(hesunit)
        END IF
     END IF   

!    ** Print hessian **


!    ** Setting the present positions & gradients as old **

     CALL DCOPY(ndf,x,1,xold,1)
     CALL DCOPY(ndf,g,1,gold,1)

!    ** Copying hessian hes to (ndf x ndf) matrix hes_mat for diagonalization **

     DO indf = 1, ndf
       DO j = 1, indf
         i = indf*(indf-1)/2 + j
         hes_mat(indf,j) = hes(i)
         hes_mat(j,indf) = hes(i)
       END DO
     END DO

     CALL diamat(hes_mat,eigvec,eigval)
  
     CALL set_hes_eig(ndf,eigval,work)

     CALL DCOPY(ndf,eigval,1,dx,1)

     CALL rat_fun_opt(ndf,dg,eigval,work,eigvec,g)

     CALL geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)

     CALL trust_radius(ndf,nfree,step,rad,rat,dr,output_unit,ionode)

      
!   ** Updating the present positions ** 

    indf = 0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list)
        natom = SIZE(atom_list)
        DO iatom=1,natom
            atom = atom_list(iatom)
            IF( gopt_param%to_move(particle_set(atom)%iatom) ) THEN
              DO j = 1, 3
                indf = indf + 1
                particle_set(atom)%r(j) = particle_set(atom)%r(j) + dr(indf)
                x(indf) = x(indf) + dr(indf)
              END DO
            END IF
        END DO
    END DO     

!  Printing section

   IF (ionode) THEN
     WRITE (UNIT=title,FMT="(A,I8,A,F20.10)") " i =",it,", E =",etot
     CALL write_particle_coordinates(particle_set,cell,dataunit,&
                                     "XMOL","POS",TRIM(title))
     CALL m_flush(dataunit)
     CALL write_qs_particle_coordinates(particle_set,cell,globenv)
     CALL write_particle_distances(particle_set,cell,globenv)
   END IF

   CALL energy_predict(ndf,work,hes,dr,g,ionode,output_unit,conv,pred)

   IF(conv)EXIT

   eold = etot

!    ** Energy & Gradients at new step **
      
     CALL qs_forces(qs_env,globenv)

     indf = 0
     DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list)
        natom = SIZE(atom_list)
        DO iatom=1,natom
            atom = atom_list(iatom)
            IF( gopt_param%to_move(particle_set(atom)%iatom) ) THEN
              DO j = 1, 3
                indf = indf + 1
                g(indf) = - particle_set(atom)%f(j)
              END DO
            END IF
        END DO
     END DO

    CALL get_qs_env(qs_env=qs_env,energy=energy)

    etot = energy%total

    ediff = etot - eold


    IF(ionode) CALL write_cycle_infos &
               (output_unit,it,etot,ediff,pred,rat,step,rad,emin)

    IF(etot < emin ) emin = etot 

    CALL update_trust_rad(rat,rad,step,ediff)

    CALL check_converg(ndf,dr,g,ionode,output_unit,conv)


    IF(conv) EXIT

  END DO

    DEALLOCATE(x)
    DEALLOCATE(hes)
    DEALLOCATE(xold)
    DEALLOCATE(g)
    DEALLOCATE(gold)
    DEALLOCATE(dx)
    DEALLOCATE(hes_mat)
    DEALLOCATE(dg)
    DEALLOCATE(eigval)
    DEALLOCATE(eigvec)
    DEALLOCATE(work)
    DEALLOCATE(dr)
    DEALLOCATE(gopt_param%to_move)

    IF(ionode)THEN
        CALL close_file(unit_number=hesunit)
        CALL close_file(unit_number=dataunit)
    END IF

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

    END SUBROUTINE geoopt_bfgs 

! *********************************************************************************

  SUBROUTINE rat_fun_opt(ndf,dg,eigval,work,eigvec,g)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: dg(ndf), eigval(ndf), &
                                                work(ndf), eigvec(ndf,ndf), &
                                                g(ndf)

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp

    INTEGER                                  :: indf, iref, iter, maxit
    LOGICAL                                  :: bisec, fail, set
    REAL(KIND=dp)                            :: DDOT, fun, fun1, fun2, fun3, &
                                                fung, lam1, lam2, ln, lp, &
                                                ssize, step, stol

!   ** Locals **

    stol  = 1.0E-8_dp
    ssize = 0.2_dp
    maxit = 999
    fail  = .FALSE.
    bisec = .FALSE.

    DO indf = 1, ndf
      dg(indf) =  DDOT(ndf,eigvec(1,indf),1,g,1)
    END DO

    set = .FALSE.

    DO 

!   calculating Lamda 

      lp = 0.0_dp
      iref = 1
      ln = 0.0_dp
      IF(eigval(iref) < 0.0_dp )ln = eigval(iref) - 0.01_dp

      iter = 0
      DO 
         iter = iter + 1
         fun = 0.0_dp
         fung = 0.0_dp
         DO indf = 1, ndf
            fun = fun + dg(indf)**2/(ln-eigval(indf))
            fung = fung - dg(indf)**2 /(ln-eigval(indf)**2)
         END DO
         fun = fun - ln
         fung = fung - 1.0_dp
         step = fun/fung
         ln = ln - step
         IF(ABS(step) < stol )GOTO 200
         IF(iter >= maxit)EXIT
      END DO
100   CONTINUE
      bisec = .TRUE.
      iter = 0
      maxit = 9999
      lam1 = 0.0_dp
      IF(eigval(iref) < 0.0_dp )lam1 = eigval(iref) - 0.01_dp
      fun1 = 0.0_dp
      DO indf = 1, ndf
        fun1 = fun1 + dg(indf)**2 / (lam1-eigval(indf))
      END DO
      fun1 = fun1 - lam1
      step = ABS(lam1)/1000.0_dp
      IF(step < ssize) step = ssize
      DO
         iter = iter + 1
         IF(iter > maxit)THEN
           ln = 0.0_dp
           lp = 0.0_dp
           fail = .TRUE.
           GOTO 300
         END IF 
         fun2 = 0.0_dp
         lam2 = lam1 - iter * step
         DO indf = 1, ndf
           fun2 = fun2 + eigval(indf)**2 / (lam2 - eigval(indf))
         END DO
         fun2 = fun2 - lam2
         IF(fun2*fun1 < 0.0_dp )THEN
           iter = 0
           DO 
             iter = iter + 1
             IF(iter > maxit)THEN
               ln = 0.0_dp
               lp = 0.0_dp
               fail = .TRUE.
               GOTO 300
             END IF
             step = (lam1+lam2)/2
             fun3 = 0.0_dp
             DO indf = 1, ndf
              fun3 = fun3 + dg(indf)**2/(step-eigval(indf))
             END DO
             fun3 = fun3 - step
           
             IF( ABS(step-lam2) < stol)THEN
               ln = step
               GOTO 200
             END IF

             IF( fun3 * fun1 < stol )THEN
                lam2 = step
             ELSE
                lam1 = step
             END IF
           END DO 
         END IF
     END DO

200  CONTINUE
     IF( (ln > eigval(iref)).OR. ( (ln > 0.0_dp) .AND. &
              (eigval(iref) > 0.0_dp) ) )THEN

         IF(.NOT. bisec)GOTO 100
         ln = 0.0_dp
         lp = 0.0_dp
         fail = .TRUE.
     END IF

300  CONTINUE

     IF(fail .AND. .NOT. set )THEN
        set = .TRUE.
        DO indf = 1, ndf
           eigval(indf) = eigval(indf) * work(indf)
        END DO
        GOTO 400
     END IF

     IF(.NOT. set) THEN
        ! CALL DCOPY(ndf,1.0,0,work,1) ! bug
        work(1:ndf)=1.0_dp
     ENDIF
      
     DO indf = 1, ndf
       eigval(indf) = eigval(indf) - ln
     END DO
     EXIT
400  END DO

  END SUBROUTINE  rat_fun_opt

! ************************************************************************************

  SUBROUTINE get_gopt_atmfix(globenv,natom,nfree,to_move)

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(IN)                      :: natom
    INTEGER, INTENT(OUT)                     :: nfree
    LOGICAL, INTENT(INOUT)                   :: to_move(natom)

    CHARACTER(len=*), PARAMETER :: routineN = 'get_gopt_atmfix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword
    INTEGER                                  :: i, itemp, ntemp
    LOGICAL                                  :: failure, section_found

  failure=.FALSE.

  CALL start_parser(file_name=globenv%input_file_name,&
       globenv=globenv,&
       start_section_label="ATMFIX",&
       section_found=section_found)

  IF (section_found) THEN
     DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
           CALL read_object(keyword)
           CALL stop_parser(routine=routineP,message="INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
           EXIT
        ELSE IF (test_result == "EOF") THEN
           CALL stop_parser(routine=routineP,message="EOF")
        ELSE
           CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
           SELECT CASE (TRIM(keyword))
           CASE ('NFIXED')
              CALL read_object(itemp)
              nfree = natom - itemp 
              DO i = 1, itemp
                 CALL read_object(ntemp,newline=.TRUE.)
                 to_move(ntemp) = .FALSE.
              END DO 
           CASE DEFAULT
              CALL stop_parser(routine=routineP,message="INVALID_KEYWORD")
           END SELECT
        END IF
    END DO
  END IF

  CALL finish_parser()

  END SUBROUTINE get_gopt_atmfix

! ******************************************************************************
 
  SUBROUTINE bfgs(ndf,dx,dg,hes,work)

!   ** Globals **

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: dx(ndf), dg(ndf), &
                                                hes(ndf*(ndf+1)/2), work(ndf)

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i, indf, j
    REAL(KIND=dp)                            :: DDOT, dxw, gdx

    CALL DSPMV('U',ndf,one,hes,dx,1,zero,work,1)

    gdx = DDOT(ndf,dg,1,dx,1)
    gdx = one/gdx
    dxw = DDOT(ndf,dx,1,work,1)
    dxw = one/dxw

    DO indf = 1, ndf
      DO j = 1, indf
         i = indf * (indf-1)/2 + j
         hes(i) = hes(i) + gdx * dg(indf) * dg(j) - &
                           dxw * work(indf) * work(j)
      END DO
    END DO

  END SUBROUTINE bfgs

!  *******************************************************************************

  SUBROUTINE set_hes_eig(ndf,eigval,work)

!    ** Globals **
    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: eigval(ndf), work(ndf)

    REAL(KIND=dp), PARAMETER                 :: max_neg = -0.5_dp, &
                                                max_pos = 5.0_dp, &
                                                min_eig = 0.005_dp, &
                                                one = 1.0_dp

    INTEGER                                  :: indf
    LOGICAL                                  :: neg

!    ** Locals **

     DO indf = 1, ndf
       IF(eigval(indf) < 0.0_dp )neg = .TRUE.
       IF(eigval(indf) > 1000.0_dp)eigval(indf) = 1000.0_dp
    END DO
    DO indf = 1, ndf
       IF( eigval(indf) < 0.0_dp ) THEN
          IF(eigval(indf) < max_neg)THEN
             eigval(indf) = max_neg
          ELSE IF(eigval(indf) > - min_eig )THEN
             eigval(indf) = - min_eig
          END IF
       ELSE IF(eigval(indf) < 1000.0_dp)THEN
          IF(eigval(indf) < min_eig)THEN
             eigval(indf) = min_eig
          ELSE IF(eigval(indf) > max_pos)THEN
             eigval(indf) = max_pos
          END IF
       END IF
     END DO

     DO indf = 1, ndf
       IF( eigval(indf) < 0.0_dp )THEN
         work(indf) = - one
       ELSE
         work(indf) =  one
       END IF
     END DO

  END SUBROUTINE set_hes_eig

! ******************************************************************************

  SUBROUTINE geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: eigval(ndf), eigvec(ndf,ndf), &
                                                hes_mat(ndf,ndf), dr(ndf), &
                                                g(ndf)

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i, indf, j
    REAL(KIND=dp), ALLOCATABLE               :: tmp(:,:)

     DO indf = 1, ndf
       eigval(indf) = one / eigval(indf)
       DO j = 1, ndf
          hes_mat(j,indf) = eigvec(j,indf) * eigval(indf)
       END DO
     END DO

     ALLOCATE(tmp(ndf,ndf))
     tmp = TRANSPOSE(eigvec)
     DO i = 1, ndf
       DO j = 1, ndf
         eigvec(i,j) = tmp(i,j)
       END DO
     END DO
     tmp = MATMUL(hes_mat,eigvec)
     DO i = 1, ndf
       DO j = 1, ndf
          hes_mat(i,j) = tmp(i,j)
       END DO
     END DO
     DEALLOCATE(tmp)

     DO i = 1, ndf
       DO j = 1, ndf
         hes_mat(j,i) = 0.5_dp * (hes_mat(j,i)+hes_mat(i,j))
         hes_mat(i,j) = hes_mat(j,i)
       END DO
     END DO

!    ** New step **

     CALL DGEMV('N',ndf,ndf,-one,hes_mat,ndf,g,1,zero,dr,1)

  END SUBROUTINE geoopt_get_step

! *****************************************************************************

  SUBROUTINE trust_radius(ndf,nfree,step,rad,rat,dr,output_unit,ionode)

!   ** Globals **
    INTEGER, INTENT(IN)                      :: ndf, nfree
    REAL(KIND=dp), INTENT(INOUT)             :: step, rad, rat, dr(ndf)
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: ionode

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp 

    REAL(KIND=dp)                            :: DDOT, scal

!   ** Locals **

    step = DDOT(ndf,dr,1,dr,1)
    scal = MAX(one,LOG10(REAL(nfree,dp)))
    step = SQRT(step/scal)

    IF( step > rad )THEN
       rat = rad / step
       CALL DSCAL(ndf,rat,dr,1)
       step = DDOT(ndf,dr,1,dr,1)
       step = SQRT(step/scal)
       IF(ionode) &
       WRITE(unit=output_unit,FMT="(/,T2,A,F8.5)") &
           " Step is scaled; Scaling factor = ", rat
       CALL m_flush(output_unit)
    END IF

  END SUBROUTINE trust_radius

! ******************************************************************************

   SUBROUTINE energy_predict(ndf,work,hes,dr,g,ionode,output_unit,conv,pred)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: work(ndf), &
                                                hes(ndf*(ndf+1)/2), dr(ndf), &
                                                g(ndf)
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(INOUT)                   :: conv
    REAL(KIND=dp), INTENT(INOUT)             :: pred

    REAL(KIND=dp), PARAMETER                 :: limit = 1.0E-9_dp , &
                                                one = 1.0_dp, zero = 0.0_dp

    REAL(KIND=dp)                            :: DDOT, ener1, ener2

     ener1 = DDOT(ndf,g,1,dr,1)
     CALL DSPMV('U',ndf,one,hes,dr,1,zero,work,1)
     ener2 = DDOT(ndf,dr,1,work,1)
     pred  = ener1 + 0.5_dp * ener2

!    ** Conv. in energy **

     conv = .FALSE.
     IF( ((pred < 0.0_dp) .AND. (ABS(pred) < limit)) .OR. &
         (ABS(pred) < 4.0E-10_dp) )THEN
        IF(ionode)THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A,F15.10)") &
            " Predicted energy difference is too small: ", pred
          WRITE(UNIT=output_unit,FMT="(/,T2,80('*'))")
          WRITE(UNIT=output_unit,FMT="(T2,2('*'),21X,A,21X,2('*'))")&
             "GEOMETRY OPTIMIZATION COMPLETED"
          WRITE(UNIT=output_unit,FMT="(T2,80('*'))")
        END  IF
        conv = .TRUE.
   END IF

  END SUBROUTINE energy_predict

! ***********************************************************************************

 SUBROUTINE write_cycle_infos(output_unit,it,etot,ediff,pred,rat,step,rad,emin)

    INTEGER, INTENT(IN)                      :: output_unit, it
    REAL(KIND=dp), INTENT(IN)                :: etot, ediff, pred, rat, step, &
                                                rad, emin

     WRITE(UNIT=output_unit,FMT="(/,T2,8('-'),A,I5,1X,8('-'))")&
     "  Informations at step = ",it
     WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
      " Total Energy               =",etot
     WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
      " Real energy change         =",ediff
     WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
      " Predicted change in energy = ",pred
     WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
      " Scaling factor             = ",rat
     WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
      " Step size                  = ",step
     WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
      " Trust radius               = ",rad
     IF( etot < emin )THEN
      WRITE(UNIT=output_unit,FMT="(T2,2A)")&
       " Decrease in energy         = ",&
       "             YES"
     ELSE
      WRITE(UNIT=output_unit,FMT="(T2,2A)")&
       " Decrease in energy         = ",&
       "              NO"
     END IF

 END SUBROUTINE write_cycle_infos

! ***********************************************************************************

   SUBROUTINE update_trust_rad(rat,rad,step,ediff)

    REAL(KIND=dp), INTENT(INOUT)             :: rat, rad, step, ediff

    REAL(KIND=dp), PARAMETER                 :: max_trust = 1.0_dp, &
                                                min_trust = 0.1_dp

    IF( rat > 4.0_dp )THEN
      IF( ediff < 0.0_dp)THEN
        rad = step * 0.5_dp
      ELSE
        rad = step * 0.25_dp
      END IF
    ELSE IF ( rat > 2.0_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 0.75_dp
      ELSE
        rad = step * 0.5_dp
      END IF
    ELSE IF ( rat > 4.0_dp/3.0_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step
      ELSE
        rad = step * 0.75_dp
      END IF
    ELSE IF (rat > 10.0_dp/9.0_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 1.25_dp
      ELSE
        rad = step
      END IF
    ELSE IF( rat > 0.9_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 1.5_dp
      ELSE
        rad = step * 1.25_dp
      END IF
    ELSE IF( rat > 0.75_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 1.25_dp
      ELSE
        rad = step
      END IF
    ELSE IF( rat > 0.5_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step
      ELSE
        rad = step * 0.75_dp
      END IF
    ELSE IF( rat > 0.25_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 0.75_dp
      ELSE
        rad = step * 0.5_dp
      END IF
    ELSE IF( ediff < 0.0_dp )THEN
        rad = step * 0.5_dp
    ELSE
        rad = step * 0.25_dp
    END IF

    rad = MAX(rad,min_trust)
    rad = MIN(rad,max_trust)

  END SUBROUTINE update_trust_rad

! ***************************************************************************

  SUBROUTINE check_converg(ndf,dr,g,ionode,output_unit,conv)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(IN)                :: dr(ndf), g(ndf)
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(OUT)                     :: conv

    REAL(KIND=dp), PARAMETER                 :: dxcon = 0.0030_dp, &
                                                gcon = 0.00045_dp, &
                                                rmsgcon = 0.00030_dp, &
                                                rmsxcon = 0.0015_dp

    INTEGER                                  :: indf
    LOGICAL                                  :: conv_dx, conv_g, conv_rdx, &
                                                conv_rg
    REAL(KIND=dp)                            :: dumm, maxdum(4)

   conv     = .FALSE.
   conv_dx  = .TRUE.
   conv_rdx = .TRUE.
   conv_g   = .TRUE.
   conv_rg  = .TRUE.

   dumm     = 0.0_dp
   DO indf = 1, ndf
      IF(indf == 1)maxdum(1) = dr(indf)
      dumm = dumm + dr(indf)**2
      IF(ABS(dr(indf)) > dxcon )conv_dx = .FALSE.
      IF(ABS(dr(indf)) > maxdum(1) )maxdum(1) = ABS(dr(indf))
   END DO
   IF( SQRT(dumm/ndf) > rmsxcon )conv_rdx = .FALSE.
   maxdum(2) = SQRT(dumm/ndf)

   dumm = 0.0_dp
   DO indf = 1, ndf
      IF(indf == 1)maxdum(3) = dr(indf)
      dumm = dumm + g(indf)**2
      IF(ABS(g(indf)) > gcon )conv_g = .FALSE.
      IF(ABS(g(indf)) > maxdum(3) )maxdum(3) = ABS(g(indf))
   END DO
   IF( SQRT(dumm/ndf) > rmsgcon )conv_rg = .FALSE.
   maxdum(4) = SQRT(dumm/ndf)

    IF(ionode)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
          " Convergence check :"
          WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Max. step size             = ",maxdum(1)
          WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Conv. limit for step size  = ",dxcon
       IF(conv_dx)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in step size   = ",&
          "             YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in step size   = ",&
          "              NO"
       END IF
          WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " RMS step size              = ",maxdum(2)
          WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Conv. limit for RMS step   = ",rmsxcon
       IF(conv_rdx)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in RMS step    = ",&
          "             YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in RMS step    = ",&
          "              NO"
       END IF
         WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Max. gradient              = ",maxdum(3)
          WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Conv. limit for gradients  = ",gcon
       IF(conv_g)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. in gradients         = ",&
          "             YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. for gradients        = ",&
          "              NO"
       END IF
         WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Max. RMS gradient          = ",maxdum(4)
          WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
          " Conv. limit for RMS grad.  = ",rmsgcon
       IF(conv_rg)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. in RMS gradients     = ",&
          "             YES"

       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. for gradients        = ",&
          "              NO"
       END IF
       WRITE(UNIT=output_unit,FMT="(T2,47('-'))")

    END IF

    IF(conv_dx .AND. conv_rdx .AND. conv_g .AND. conv_rg)conv = .TRUE.

    IF(conv)THEN
      IF(ionode)THEN
        WRITE(UNIT=output_unit,FMT="(/,T2,80('*'))")
        WRITE(UNIT=output_unit,FMT="(T2,2('*'),21X,A,21X,2('*'))")&
           "GEOMETRY OPTIMIZATION COMPLETED"
        WRITE(UNIT=output_unit,FMT="(T2,80('*'))")
      END IF
    END IF

  END SUBROUTINE check_converg

END MODULE bfgs_optimizer 
