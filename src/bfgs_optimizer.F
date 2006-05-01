!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/bfgs_optimizer [1.0] *
!!
!!   NAME
!!     bfgs_optimizer
!!
!!   FUNCTION
!!     Routines for Geometry optimization using BFGS algorithm 
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************

MODULE bfgs_optimizer

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_param_type,&
                                             gopt_f_type
  USE input_constants,                 ONLY: dump_xmol
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: geoopt_bfgs,&
            write_cycle_infos,&
            check_converg

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bfgs_optimizer'


CONTAINS

!!*** **********************************************************************

  SUBROUTINE geoopt_bfgs(force_env,gopt_param,globenv,force_env_section,&
       gopt_env, x0, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'geoopt_bfgs'
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    CHARACTER(LEN=80)                        :: title
    INTEGER :: atom, dataunit, handle, hesunit, i, iatom, ihandle, ikind, &
      indf, istat, it, j, maxiter, natom, ndf, nfixed, nfree, nkind, &
      output_unit
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: conv, hesrest, ionode, &
                                                should_stop
    REAL(KIND=dp)                            :: ediff, emin, eold, etot, &
                                                pred, rad, rat, step
    REAL(KIND=dp), ALLOCATABLE :: dg(:), dr(:), dx(:), eigval(:), &
      eigvec(:,:), g(:), gold(:), hes(:), hes_mat(:,:), work(:), x(:), xold(:)
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL cp_add_iter_level(logger%iter_info,"GEO_OPT",ihandle,error)
    CALL timeset(routineN,"I","",handle)

    ionode      = globenv%ionode
    maxiter     = gopt_param%max_iter
    hesrest     = gopt_param%hesrest
    conv        = .FALSE.
    rad         = 0.10_dp
    rat         = 0.0_dp

    output_unit = cp_print_key_unit_nr(logger,globenv%input_file,&
         "MOTION%GEOOPT%PRINT%PROGRAM_RUN_INFO",&
         extension=".geoLog",error=error)

    CALL force_env_get(force_env,&
         cell=cell,&
         subsys=subsys)

    CALL cp_subsys_get(subsys=subsys(1)%subsys,&
         atomic_kinds=atomic_kinds,&
         particles=particles)

    nkind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els

    natom = particles%n_els  
    particle_set => particles%els    

    !   Atom fixing constraints for optimizer

    nfixed = COUNT(particle_set(:)%is_fixed)
    nfree = natom - nfixed
    ndf = 3*nfree

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
          IF (.NOT.particle_set(atom)%is_fixed) THEN
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

    IF(output_unit>0)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       WRITE(UNIT=output_unit,FMT="(T2,A,T28,A,T78,A)")&
            "***","STARTING GEOMETRY OPTIMIZATION","***"
       WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       WRITE(UNIT=output_unit,FMT="(/,(T2,A,I10))")&
            "Number of free atoms         = ",nfree,&
            "Number of degrees of freedom = ",ndf
       CALL m_flush(output_unit)
    END IF

    CALL write_qs_particle_coordinates(particle_set,cell,force_env_section,error) 
    CALL write_particle_distances(particle_set,cell,force_env_section,error)
    CALL write_structure_data(particle_set,cell,force_env_section,error)

    IF(ionode)THEN
       CALL open_file(file_name="HESSIAN",&
            file_status="UNKNOWN",&
            file_form="UNFORMATTED",&
            file_action="READWRITE",&
            unit_number=hesunit)
    END IF

    ! Calculate Energy & Gradients

    CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,error=error)

    ! Copy forces to g, the gradients

    indf = 0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            atom_list=atom_list)
       natom = SIZE(atom_list)
       DO iatom=1,natom
          atom = atom_list(iatom)
          IF (.NOT.particle_set(atom)%is_fixed) THEN
             DO j = 1, 3
                indf = indf + 1
                g(indf) = -particle_set(atom)%f(j)
             END DO
          END IF
       END DO
    END DO

    CALL force_env_get(force_env=force_env,potential_energy=etot)

    emin = etot

    it = 0

    DO 
       
       it = it + 1
       CALL cp_iterate(logger%iter_info,iter_nr=it,error=error)
       IF(it > maxiter )THEN
          IF(output_unit>0)THEN
             WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                  "*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***"
             WRITE(UNIT=output_unit,FMT="(T2,A)")&
                  "***        EXITING GEOMETRY OPTIMIZATION         ***" 
             CALL m_flush(output_unit)
          END IF
          EXIT
       END IF
       ! check for an external exit command 
       CALL external_control(should_stop,"GEO",globenv,error)
       IF(should_stop) EXIT
       CALL write_restart(force_env=force_env,globenv=globenv,force_env_section=force_env_section)

       IF (output_unit>0)THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,26('-'))")
          WRITE(UNIT=output_unit,FMT="(T2,A,I6)")&
               "OPTIMIZATION STEP: ",it
          WRITE(UNIT=output_unit,FMT="(T2,26('-'))")
          CALL m_flush(output_unit)
       END IF

       !    ** Hessian update/restarting **

       IF ((it == 1).AND.hesrest) THEN
          IF(ionode)THEN
             DO i=1,SIZE(hes)
                READ(UNIT=hesunit,IOSTAT=istat)hes(i)
             END DO
             REWIND(hesunit)
             IF (istat /= 0) THEN
                WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                     "**  Error while reading HESSIAN **"
                conv = .TRUE.
             ELSE
                WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                     "*** Initial Hessian has been restarted successfully ***"
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

       !    ** Setting the present positions & gradients as old **

       CALL DCOPY(ndf,x,1,xold,1)
       CALL DCOPY(ndf,g,1,gold,1)

       !    ** Copying hessian hes to (ndf x ndf) matrix hes_mat for diagonalization **

       DO indf = 1, ndf
          DO j = 1, indf
             i = indf*(indf-1)/2 + j
             eigvec(indf,j) = hes(i)
             eigvec(j,indf) = hes(i)
          END DO
       END DO

       CALL diamat_all(eigvec,eigval,dac=.FALSE.)

       CALL set_hes_eig(ndf,eigval,work)

       CALL DCOPY(ndf,eigval,1,dx,1)

       CALL rat_fun_opt(ndf,dg,eigval,work,eigvec,g)

       CALL geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)

       CALL trust_radius(ndf,nfree,step,rad,rat,dr,output_unit,output_unit>0)

       !   ** Updating the present positions ** 

       indf = 0
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               atom_list=atom_list)
          natom = SIZE(atom_list)
          DO iatom=1,natom
             atom = atom_list(iatom)
             IF (.NOT.particle_set(atom)%is_fixed) THEN
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
          dataunit = cp_print_key_unit_nr(logger,globenv%input_file,&
               "MOTION%PRINT%TRAJECTORY",extension=".xyz",error=error)
          IF (dataunit>0) THEN
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)") " i =",it,", E =",etot
             CALL write_particle_coordinates(particle_set,cell,dataunit,&
                  dump_xmol,"POS",TRIM(title))
             CALL m_flush(dataunit)
          END IF
          CALL cp_print_key_finished_output(dataunit,logger,globenv%input_file,&
               "MOTION%PRINT%TRAJECTORY", error=error)     
          CALL write_qs_particle_coordinates(particle_set,cell,force_env_section,error)
       END IF

       CALL write_particle_distances(particle_set,cell,force_env_section,error)
       CALL write_structure_data(particle_set,cell,force_env_section,error)

       CALL energy_predict(ndf,work,hes,dr,g,output_unit>0,output_unit,conv,pred)

       IF (conv) EXIT

       eold = etot

       !    ** Energy & Gradients at new step **

       CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,error=error)

       indf = 0
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               atom_list=atom_list)
          natom = SIZE(atom_list)
          DO iatom=1,natom
             atom = atom_list(iatom)
             IF (.NOT.particle_set(atom)%is_fixed) THEN
                DO j = 1, 3
                   indf = indf + 1
                   g(indf) = -particle_set(atom)%f(j)
                END DO
             END IF
          END DO
       END DO

       CALL force_env_get(force_env=force_env,potential_energy=etot)

       ediff = etot - eold

       IF(output_unit>0) CALL write_cycle_infos &
            (output_unit,it,etot,ediff,pred,rat,step,rad,emin," BFGS")

       IF(etot < emin ) emin = etot 

       CALL update_trust_rad(rat,rad,step,ediff)

       CALL check_converg(ndf,dr,g,output_unit>0,output_unit,conv,gopt_param)


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

    IF(ionode)THEN
       CALL close_file(unit_number=hesunit)
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,globenv%input_file,&
         "MOTION%GEOOPT%PRINT%PROGRAM_RUN_INFO", error=error)     

    CALL cp_rm_iter_level(logger%iter_info,ihandle,error)
    CALL timestop(0.0_dp,handle)

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
    REAL(KIND=dp)                            :: tmp(ndf,ndf)

     hes_mat(:,:) = eigvec(:,:)

     DO indf=1,ndf
       eigval(indf) = one/eigval(indf)
       CALL dscal(ndf,eigval(indf),hes_mat(1,indf),1)
     END DO

     CALL dgemm ("N","T",ndf,ndf,ndf,one,hes_mat,ndf,eigvec,ndf,zero,&
                 tmp,ndf)

     DO i=1,ndf
       DO j=1,ndf
         hes_mat(j,i) = 0.5_dp*(tmp(j,i) + tmp(i,j))
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
       IF(ionode) THEN
         WRITE(unit=output_unit,FMT="(/,T2,A,F8.5)") &
             " Step is scaled; Scaling factor = ", rat
         CALL m_flush(output_unit)
       ENDIF
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

 SUBROUTINE write_cycle_infos(output_unit,it,etot,ediff,pred,rat,step,rad,emin,wildcard)

    INTEGER, INTENT(IN)                      :: output_unit, it
    REAL(KIND=dp), INTENT(IN)                :: etot, ediff
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: pred, rat, step, rad
    REAL(KIND=dp), INTENT(IN)                :: emin
    CHARACTER(LEN=5), OPTIONAL               :: wildcard

     WRITE(UNIT=output_unit,FMT="(/,T2,8('-'),A,I5,1X,12('-'))")&
     "  Informations at step = ",it
     WRITE(UNIT=output_unit,FMT="(T2,A,T47,A)")&
      " Optimization Method        = ",wildcard
     WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
      " Total Energy               = ",etot
     WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
      " Real energy change         = ",ediff
     IF (PRESENT(pred)) THEN
        WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
             " Predicted change in energy = ",pred
     END IF
     IF (PRESENT(rat)) THEN
        WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
             " Scaling factor             = ",rat
     END IF
     IF (PRESENT(step)) THEN
        WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
             " Step size                  = ",step
     END IF
     IF (PRESENT(rad)) THEN
        WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
             " Trust radius               = ",rad
     END IF
     IF( etot < emin )THEN
      WRITE(UNIT=output_unit,FMT="(T2,2A)")&
       " Decrease in energy         =     ",&
       "             YES"
     ELSE
      WRITE(UNIT=output_unit,FMT="(T2,2A)")&
       " Decrease in energy         =     ",&
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

  SUBROUTINE check_converg(ndf,dr,g,ionode,output_unit,conv,gopt_param)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(IN)                :: dr(ndf), g(ndf)
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(OUT)                     :: conv
    TYPE(gopt_param_type), POINTER           :: gopt_param

    INTEGER                                  :: indf
    LOGICAL                                  :: conv_dx, conv_g, conv_rdx, &
                                                conv_rg
    REAL(KIND=dp)                            :: dumm, dxcon, gcon, maxdum(4), &
                                                rmsgcon, rmsxcon

!    REAL(KIND=dp), PARAMETER                 :: dxcon = 0.0030_dp, &
!                                                gcon = 0.00045_dp, &
!                                                rmsgcon = 0.00030_dp, &
!                                                rmsxcon = 0.0015_dp

    dxcon = gopt_param%max_dr
    gcon = gopt_param%max_force
    rmsgcon = gopt_param%rms_force
    rmsxcon = gopt_param%rms_dr

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
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Max. step size             = ",maxdum(1)
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Conv. limit for step size  = ",dxcon
       IF(conv_dx)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in step size   =     ",&
          "             YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in step size   =     ",&
          "              NO"
       END IF
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " RMS step size              = ",maxdum(2)
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Conv. limit for RMS step   = ",rmsxcon
       IF(conv_rdx)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in RMS step    =     ",&
          "             YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Convergence in RMS step    =     ",&
          "              NO"
       END IF
         WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Max. gradient              = ",maxdum(3)
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Conv. limit for gradients  = ",gcon
       IF(conv_g)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. in gradients         =     ",&
          "             YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. for gradients        =     ",&
          "              NO"
       END IF
         WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Max. RMS gradient          = ",maxdum(4)
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
          " Conv. limit for RMS grad.  = ",rmsgcon
       IF(conv_rg)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. in RMS gradients     =     ",&
          "             YES"

       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
          " Conv. for gradients        =     ",&
          "              NO"
       END IF
       WRITE(UNIT=output_unit,FMT="(T2,51('-'))")

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
! *****************************************************************************


END MODULE bfgs_optimizer 
