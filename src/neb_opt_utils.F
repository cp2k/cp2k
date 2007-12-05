!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/neb_opt_utils [1.0] *
!!
!!   NAME
!!     neb_opt_utils
!!
!!   FUNCTION
!!     Module with utility to perform MD Nudged Elastic Band Calculation
!!
!!   NOTES
!!     Numerical accuracy for parallel runs:
!!      Each replica starts the SCF run from the one optimized
!!      in a previous run. It may happen then energies and derivatives
!!      of a serial run and a parallel run could be slightly different
!!      'cause of a different starting density matrix.
!!      Exact results are obtained using:
!!         EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!!
!!   AUTHOR
!!     Teodoro Laino 10.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
MODULE neb_opt_utils
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE neb_types,                       ONLY: neb_type
  USE neb_utils,                       ONLY: neb_calc_energy_forces,&
                                             reorient_images
  USE particle_types,                  ONLY: particle_type
  USE replica_types,                   ONLY: replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_opt_utils'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

  PUBLIC :: accept_diis_step,&
            neb_ls
            
  REAL(KIND=dp), DIMENSION(2:10), PRIVATE :: acceptance_factor = &
       (/0.97_dp, 0.84_dp, 0.71_dp, 0.67_dp, 0.62_dp, 0.56_dp, 0.49_dp, 0.41_dp, 0.0_dp/)

CONTAINS

!!****f* neb_opt_utils/accept_diis_step [1.0] *
!!
!!   NAME
!!     accept_diis_step
!!
!!   FUNCTION
!!     Performes few basic operations for the NEB DIIS optimization
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 10.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  FUNCTION accept_diis_step(apply_diis, n_diis,err,crr,set_err,sline,coords,&
       output_unit,check_diis,error) RESULT(accepted)
    LOGICAL, INTENT(IN)                      :: apply_diis
    INTEGER, INTENT(IN)                      :: n_diis
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: err, crr
    INTEGER, DIMENSION(:), POINTER           :: set_err
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sline, coords
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: check_diis
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: accepted

    CHARACTER(len=*), PARAMETER :: routineN = 'accept_diis_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iend, ind, indi, indj, &
                                                info, istart, iv, iw, j, jv, &
                                                k, lwork, np, nv, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: IWORK
    LOGICAL                                  :: failure, increase_error, &
                                                ionode
    REAL(dp), DIMENSION(:, :), POINTER       :: work
    REAL(KIND=dp)                            :: eps_svd
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: S, Work_dgesdd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: U, VT, wrk, wrk_inv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: awrk, cwrk, ref, step
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode.AND.debug_this_module) THEN
      iw = cp_logger_get_default_unit_nr(logger)
    ELSE
      iw = -1
    END IF
    failure  = .FALSE.
    accepted = .FALSE.
    IF (.NOT.failure) THEN
       ! find the index with the minimum element of the set_err array 
       nv = MINLOC(set_err,1)
       IF (iw>0) WRITE(iw,'(A,I3)')"Entering into the DIIS module. Error vector number:: ",nv
       set_err(nv) = 1
       eps_svd     = 1.0E-10_dp
       ALLOCATE(step(SIZE(sline,1)*SIZE(sline,2)),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ref(SIZE(sline,1)*SIZE(sline,2)),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       err(:,nv) = RESHAPE(sline,(/SIZE(sline,1)*SIZE(sline,2)/))
       crr(:,nv) = RESHAPE(coords,(/SIZE(sline,1)*SIZE(sline,2)/))
       jv = n_diis
       IF (ALL(set_err==1).AND.apply_diis) THEN
          IF (iw>0) WRITE(iw,'(A)')"Applying DIIS equations"
          ! Apply DIIS..
          DO jv = 2, n_diis
             np = jv+1
             IF (iw>0) WRITE(iw,'(A)')"Applying DIIS equations with the last"//&
                  cp_to_string(jv)//" error vectors"
             ALLOCATE(wrk(np,np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(work(np,np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(wrk_inv(np,np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(cwrk(np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(awrk(np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             awrk           = 0.0_dp
             wrk            = 1.0_dp
             wrk(np,np)     = 0.0_dp
             awrk(np)       = 1.0_dp
             DO i = 1, jv
                indi = n_diis - i + 1
                DO j = i, jv
                   indj = n_diis - j + 1
                   wrk(i,j) = DOT_PRODUCT(err(:,indi),err(:,indj))
                   wrk(j,i) = wrk(i,j)
                END DO
             END DO
             IF (iw>0) WRITE(iw,'(A)')"DIIS Matrix"//cp_to_string(np)//"x"//cp_to_string(np)//"."
             IF (iw>0) WRITE(iw,'('//cp_to_string(np)//'F12.6)')wrk
             ! Inverte the DIIS Matrix
             work=TRANSPOSE(wrk)
             ! Workspace query
             ALLOCATE(iwork(8*np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(S(np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(U(np,np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(VT(np,np),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(work_dgesdd(1),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             lwork=-1
             CALL DGESDD('S',np,np,work,np,S,U,np,vt,np,work_dgesdd,lwork,iwork,info)
             lwork=work_dgesdd(1)
             DEALLOCATE(work_dgesdd,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(work_dgesdd(lwork),stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL DGESDD('S',np,np,work,np,S,U,np,vt,np,work_dgesdd,lwork,iwork,info)
             ! Construct the inverse
             DO k=1,np
                ! Invert SV
                IF (S(k)<eps_svd) THEN
                   S(k)=0.0_dp
                ELSE
                   S(k)=1.0_dp/S(k)
                ENDIF
                VT(k,:)=VT(k,:)*S(k)
             ENDDO
             CALL BL_GEMM('T','T',np,np,np,1.0_dp,VT,np,U,np,0.0_dp,wrk_inv,np)
             DEALLOCATE(iwork,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(S,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(U,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(VT,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(work_dgesdd,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             cwrk = MATMUL(wrk_inv,awrk)
             ! Check the DIIS solution
             step = 0.0_dp
             ind  = 0
             DO i = n_diis, n_diis-jv+1, -1
                ind = ind + 1
                step = step + (crr(:,i)+err(:,i))*cwrk(ind)
             END DO
             step = step - crr(:,n_diis)
             ref  = err(:,n_diis)
             increase_error = check_diis_solution(jv,n_diis,cwrk,step,ref,&
                                                  output_unit,check_diis,error)
             ! possibly enlarge the error space
             IF (increase_error) THEN
                accepted = .TRUE.
                sline    = RESHAPE(step,(/SIZE(sline,1),SIZE(sline,2)/))
             ELSE
                DEALLOCATE(awrk,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(cwrk,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(wrk,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(work,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(wrk_inv,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)                
                EXIT
             END IF
             DEALLOCATE(awrk,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(cwrk,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(wrk,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(work,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(wrk_inv,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END DO
          IF (iw>0) WRITE(iw,'(A)')"Exiting DIIS accepting"//&
               cp_to_string(MIN(n_diis,jv))//" errors."
       END IF
       IF (ALL(set_err==1)) THEN
          ! always delete the last error vector from the history vectors
          ! move error vectors and the set_err in order to have free space
          ! at the end of the err array
          istart = MAX(2,n_diis-jv+2)
          iend   = n_diis
          indi   = 0
          DO iv = istart, iend
             indi = indi + 1
             err(:,indi)   = err(:,iv)
             crr(:,indi)   = crr(:,iv)
             set_err(indi) = 1
          END DO
          DO iv = indi+1, iend
             set_err(iv) = -1
          END DO
       END IF
    END IF
    DEALLOCATE(step,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ref,stat=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END FUNCTION accept_diis_step

!!****f* neb_opt_utils/check_diis_solution [1.0] *
!!
!!   NAME
!!     check_diis_solution
!!
!!   FUNCTION
!!     Check conditions for the acceptance of the DIIS step
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 10.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  FUNCTION check_diis_solution(nv, n_diis, cwrk, step, ref, output_unit, check_diis, error)&
   RESULT(accepted)
    INTEGER, INTENT(IN)                      :: nv, n_diis
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cwrk, step, ref
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: check_diis
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: accepted

    CHARACTER(len=*), PARAMETER :: routineN = 'check_diis_solution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: costh, norm1, norm2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tmp

    failure  = .FALSE.
    accepted = .TRUE.
    ALLOCATE(tmp(SIZE(step)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (accepted) THEN
          ! (a) The direction of the DIIS step, can be compared to the reference step.
          !     if the angle is grater than a specified value, the DIIS step is not
          !     acceptable. 
          norm1 = SQRT(DOT_PRODUCT(ref,ref))
          norm2 = SQRT(DOT_PRODUCT(step,step))
          costh = DOT_PRODUCT(ref,step)/(norm1*norm2)
          IF (check_diis) THEN
             IF (costh < acceptance_factor(MIN(10,nv))) accepted = .FALSE.
          ELSE
             IF (costh <= 0.0_dp) accepted = .FALSE.
          END IF
          IF (output_unit > 0.AND.(.NOT.accepted))  THEN
             WRITE(output_unit,'(T2,"DIIS|",A)')&
                  "The direction of the DIIS step, can be compared to the reference step.",&
                  "If the angle is grater than a specified value, the DIIS step is not",&
                  "acceptable. Value exceeded. Reset DIIS!"
             WRITE(output_unit,'(T2,"DIIS|",A,F6.3,A,F6.3,A)')&
                  "Present Cosine <",costh,"> compared with the optimal value <",&
                  acceptance_factor(MIN(10,nv)),"> ."
          END IF
       END IF
       IF (accepted.AND.check_diis) THEN
          ! (b) The length of the DIIS step is limited to be no more than 10 times
          !     the reference step
          IF (norm1 > norm2*10.0_dp) accepted = .FALSE.
          IF (output_unit > 0.AND.(.NOT.accepted))  THEN
             WRITE(output_unit,'(T2,"DIIS|",A)')&
                  "The length of the DIIS step is limited to be no more than 10 times",&
                  "the reference step. Value exceeded. Reset DIIS!"
          END IF
       END IF
       IF (accepted.AND.check_diis) THEN
          ! (d) If the DIIS matrix is nearly singular, the norm of the DIIS step
          !     vector becomes small and cwrk/norm1 becomes large, signaling a
          !     numerical stability problems. If the magnitude of cwrk/norm1 
          !     exceeds 10^8 then the step size is assumed to be unacceptable
          IF (ANY(ABS(cwrk(1:nv)/norm1)>10**8_dp)) accepted = .FALSE.
          IF (output_unit > 0.AND.(.NOT.accepted))  THEN
             WRITE(output_unit,'(T2,"DIIS|",A)')&
                  "If the DIIS matrix is nearly singular, the norm of the DIIS step",&
                  "vector becomes small and Coeff/E_norm becomes large, signaling a",&
                  "numerical stability problems. IF the magnitude of Coeff/E_norm",&
                  "exceeds 10^8 THEN the step size is assumed to be unacceptable.",&
                  "Value exceeded. Reset DIIS!"
          END IF          
       END IF
    END IF
    DEALLOCATE(tmp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END FUNCTION check_diis_solution

!!****f* neb_opt_utils/neb_ls [1.0] *
!!
!!   NAME
!!     neb_ls
!!
!!   FUNCTION
!!     Perform a line minimization search in optimizing a NEB with DIIS
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 10.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE  neb_ls(stepsize, sline, rep_env, neb_env, coords, energies, forces,&
               vels, particle_set, iw, output_unit, distances, diis_section, &
               error)
    REAL(KIND=dp), INTENT(INOUT)             :: stepsize
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sline
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), OPTIONAL, POINTER        :: neb_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: energies
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces, vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: iw, output_unit
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: distances
    TYPE(section_vals_type), POINTER         :: diis_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_ls', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, np, stat
    LOGICAL                                  :: failure
    ! replaced xc by xc_cray to work around yet another bug in pgf90 on CRAY xt3 
    REAL(KIND=dp)                            :: a, b, max_stepsize, xa, xb, &
                                                xc_cray, ya, yb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Icoord

    failure  = .FALSE.
    IF (.NOT.failure) THEN
       ALLOCATE(Icoord(SIZE(coords,1),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL section_vals_val_get(diis_section,"NP_LS", i_val=np, error=error)
       CALL section_vals_val_get(diis_section,"MAX_STEPSIZE", r_val=max_stepsize, error=error)
       Icoord = coords
       xa = 0.0_dp
       ya = SUM(sline*forces)
       xb = xa + MIN(SUM(sline*forces)*stepsize,max_stepsize)
       xc_cray = xb
       i  = 1
       DO WHILE (i<=np-1)
          i = i + 1
          coords = Icoord + xb*sline
          CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
               output_unit, distances, rep_env%para_env, error)
          neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
          CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
               particle_set, iw, error)
          yb = SUM(sline*forces)
          a  = (ya-yb)/(2.0_dp*(xa-xb))
          b  = ya - 2.0_dp*a*xa
          xc_cray = -b/(2.0_dp*a)
          IF (xc_cray>max_stepsize) THEN
             IF (output_unit>0) WRITE(output_unit,'(T2,2(A,F6.3),A)')&
                  "LS| Predicted stepsize (",xc_cray,") greater than allowed stepsize (",&
                  max_stepsize,").  Reset!"
             xc_cray = max_stepsize
             EXIT
          END IF
          ! No Extrapolation .. only interpolation
          IF ((xc_cray<=MIN(xa,xb) .OR. xc_cray>=MAX(xa,xb)).AND.(ABS(xa-xb)>1.0E-5_dp)) THEN
             IF (output_unit>0) WRITE(output_unit,'(T2,2(A,I5),A)')&
                  "LS| Increasing the number of point from ",np," to ",np+1,"."
             np = np + 1
          END IF
          !
          IF (ABS(yb)<ABS(ya)) THEN
             ya = yb
             xa = xb
          END IF
          xb = xc_cray
       END DO
       stepsize = xc_cray
       coords = Icoord
       DEALLOCATE(Icoord,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE neb_ls

END MODULE neb_opt_utils
