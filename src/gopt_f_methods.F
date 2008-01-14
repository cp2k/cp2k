!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****h* cp2k/gopt_f_methods [1.0] *
!!
!!   NAME
!!     gopt_f_methods
!!
!!   FUNCTION
!!     contains a functional that calculates the energy and its derivatives
!!     for the geometry optimizer
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     -
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE gopt_f_methods
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE dimer_methods,                   ONLY: cp_eval_at_ts
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_minimization_method_id,&
                                             default_ts_method_id,&
                                             dump_atomic,&
                                             dump_dcd,&
                                             dump_xmol
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             sp
  USE machine,                         ONLY: m_flush
  USE md_energies,                     ONLY: get_output_format
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates,&
                                             write_structure_data
  USE qmmm_util,                       ONLY: apply_qmmm_translate
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_f_methods'

  PUBLIC :: cp_eval_at, gopt_f_create_x0
  PUBLIC :: write_cycle_infos, check_converg, write_final_info, write_geo_traj

CONTAINS

!!****f* gopt_f_methods/cp_eval_at [1.0] *
!!
!!   NAME
!!     cp_eval_at
!!
!!   FUNCTION
!!     evaluete the potential energy and its gradients using an array
!!     with same dimension as the particle_set
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment
!!     x: the position where the function should be evaluated
!!     f: the function value
!!     gradient: the value of its gradient
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_eval_at(gopt_env, x, f, gradient, master, para_env, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: x
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    INTEGER, INTENT(in)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idg, idir, ip, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: f_ts
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient_ts
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    NULLIFY(particles, subsys, cell, gradient_ts)

    CALL timeset(routineN,"I","",handle)

    CALL force_env_get(gopt_env%force_env, subsys=subsys, cell=cell,error=error)
    idg=0
    CALL cp_subsys_get(subsys, particles=particles, error=error)
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          particles%els(ip)%r(idir)=x(idg)
       END DO
    END DO
    CALL write_structure_data(particles%els,cell,gopt_env%motion_section,error)

    SELECT CASE(gopt_env%type_id)
    CASE(default_minimization_method_id)
       ! Geometry Minimization
       CALL force_env_calc_energy_force(gopt_env%force_env,&
            calc_force=PRESENT(gradient),error=error)

       ! Possibly take the potential energy
       IF (PRESENT(f)) THEN
          CALL force_env_get(gopt_env%force_env,potential_energy=f,error=error)
       END IF
       
       ! Possibly take the gradients
       IF (PRESENT(gradient)) THEN
          IF (master==para_env%mepos) THEN ! we are on the master
             CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
             CALL force_env_get(gopt_env%force_env, subsys=subsys,error=error)
             idg=0
             CALL cp_subsys_get(subsys, particles=particles, error=error)
             DO ip=1,particles%n_els
                DO idir=1,3
                   idg=idg+1
                   CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
                   gradient(idg)=-particles%els(ip)%f(idir)
                END DO
             END DO
          END IF
       END IF
    CASE(default_ts_method_id)
       ! Transition State Optimization
       IF (PRESENT(gradient)) THEN
          ALLOCATE(gradient_ts(particles%n_els*3),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       ! Real calculation of energy and forces for transition state optimization
       CALL cp_eval_at_ts (gopt_env, x, f_ts, gradient_ts, calc_force=PRESENT(gradient), error=error)

       ! Possibly take the potential energy
       IF (PRESENT(f)) f = f_ts
       
       ! Possibly take the gradients
       IF (PRESENT(gradient)) THEN
          IF (master==para_env%mepos) THEN ! we are on the master
             CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
             gradient = gradient_ts
          END IF
          DEALLOCATE(gradient_ts, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END SELECT

    ! This call is necessary for QM/MM if a Translation is applied
    ! this makes the geometry optimizer consistent
    idg=0
    CALL cp_subsys_get(subsys, particles=particles,error=error)
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          particles%els(ip)%r(idir)=x(idg)
       END DO
    END DO
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_eval_at

!!****f* gopt_f_methods/gopt_f_create_x0 [1.0] *
!!
!!   NAME
!!     gopt_f_create_x0
!!
!!   FUNCTION
!!     returns the value of the parameters for the actual configuration
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment you want the info about
!!     x0: the parameter vector (is allocated by this routine)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE gopt_f_create_x0(gopt_env, x0, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_create_x0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idg, idir, ip, ndg, stat
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    NULLIFY(particles,subsys)

    CALL force_env_get(gopt_env%force_env, subsys=subsys,&
         error=error)
    ! before starting we handle the case of translating coordinates (QM/MM)
    CALL apply_qmmm_translate(gopt_env%force_env,error=error)
    CALL cp_subsys_get(subsys, particles=particles,error=error)
    ndg=particles%n_els*3

    ALLOCATE(x0(ndg),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    idg=0
    CALL cp_subsys_get(subsys, particles=particles,error=error)
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          CPInvariant(SIZE(x0)>=idg,cp_failure_level,routineP,error,failure)
          x0(idg)=particles%els(ip)%r(idir)
       END DO
    END DO

  END SUBROUTINE gopt_f_create_x0
!***************************************************************************

!!****** gopt_f_methods/write_cycle_infos [1.0] *
!!
!!   NAME
!!    write_cycle_infos
!!
!!   FUNCTION
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
 SUBROUTINE write_cycle_infos(output_unit,it,etot,ediff,pred,rat,step,rad,emin,wildcard)

    INTEGER, INTENT(IN)                      :: output_unit, it
    REAL(KIND=dp), INTENT(IN)                :: etot
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: ediff, pred, rat, step, rad, &
                                                emin
    CHARACTER(LEN=5), OPTIONAL               :: wildcard

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_cycle_infos', &
      routineP = moduleN//':'//routineN

     WRITE(UNIT=output_unit,FMT="(/,T2,8('-'),A,I5,1X,12('-'))")&
     "  Informations at step = ",it
     WRITE(UNIT=output_unit,FMT="(T2,A,T47,A)")&
      " Optimization Method        = ",wildcard
     WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
      " Total Energy               = ",etot
     IF (PRESENT(ediff)) THEN
        WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
             " Real energy change         = ",ediff
     END IF
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
     IF (PRESENT(emin)) THEN
        IF( etot < emin )THEN
           WRITE(UNIT=output_unit,FMT="(T2,2A)")&
                " Decrease in energy         =     ",&
                "             YES"
        ELSE
           WRITE(UNIT=output_unit,FMT="(T2,2A)")&
                " Decrease in energy         =     ",&
                "              NO"
        END IF
     END IF
     IF (it==0) WRITE(UNIT=output_unit,FMT="(T2,51('-'))")
 END SUBROUTINE write_cycle_infos

!!****** gopt_f_methods/check_converg [1.0] *
!!
!!   NAME
!!    check_converg  
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
  SUBROUTINE check_converg(ndf,dr,g,ionode,output_unit,conv,gopt_param)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(IN)                :: dr(ndf), g(ndf)
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(OUT)                     :: conv
    TYPE(gopt_param_type), POINTER           :: gopt_param

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_converg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: indf
    LOGICAL                                  :: conv_dx, conv_g, conv_rdx, &
                                                conv_rg
    REAL(KIND=dp)                            :: dumm, dxcon, gcon, maxdum(4), &
                                                rmsgcon, rmsxcon

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
       IF(indf == 1)maxdum(1) = ABS(dr(indf))
       dumm = dumm + dr(indf)**2
       IF(ABS(dr(indf)) > dxcon )conv_dx = .FALSE.
       IF(ABS(dr(indf)) > maxdum(1) )maxdum(1) = ABS(dr(indf))
    END DO
    IF( SQRT(dumm/ndf) > rmsxcon )conv_rdx = .FALSE.
    maxdum(2) = SQRT(dumm/ndf)

    dumm = 0.0_dp
    DO indf = 1, ndf
       IF(indf == 1)maxdum(3) = ABS(g(indf))
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
            " RMS gradient               = ",maxdum(4)
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

    IF (conv_dx .AND. conv_rdx .AND. conv_g .AND. conv_rg) conv = .TRUE.

    IF (conv) THEN
      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
        WRITE (UNIT=output_unit,FMT="(T2,A,T25,A,T78,A)")&
          "***","GEOMETRY OPTIMIZATION COMPLETED","***"
        WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
      END IF
    END IF

  END SUBROUTINE check_converg

!!****** gopt_f_methods/write_final_info [1.0] *
!!
!!   NAME
!!    write_final_info
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!    Teodoro Laino [tlaino] - University of Zurich - 11.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
  SUBROUTINE write_final_info(conv, it, gopt_env, x0, master, para_env, force_env,&
    motion_section, root_section, force_env_section, error)
    LOGICAL, INTENT(IN)                      :: conv
    INTEGER, INTENT(INOUT)                   :: it
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    INTEGER, INTENT(IN)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: motion_section, root_section, &
                                                force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_final_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: etot
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL force_env_get(force_env, cell=cell, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys=subsys, particles=particles,error=error)
    particle_set => particles%els
    IF (conv) THEN
      it = it + 1
      CALL write_structure_data(particle_set,cell,motion_section,error)
      CALL write_restart(force_env=force_env,root_section=root_section,&
                         force_env_section=force_env_section,error=error)

      ! Write out more data for the minimum
      IF (logger%para_env%mepos==logger%para_env%source) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger), fmt="(/,T20,' Reevaluating energy at the minimum ')")
      END IF
      CALL cp_eval_at(gopt_env,x0,f=etot,master=master,para_env=para_env,error=error)
      CALL write_geo_traj(force_env,  root_section, it, etot, error)
    END IF

  END SUBROUTINE write_final_info

!!****** gopt_f_methods/write_geo_traj [1.0] *
!!
!!   NAME
!!    write_geo_traj
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!    Teodoro Laino [tlaino] - University of Zurich - 11.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
  SUBROUTINE write_geo_traj(force_env,  root_section, it, etot, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: it
    REAL(KIND=dp), INTENT(IN)                :: etot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_geo_traj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_act, my_ext, my_form, &
                                                my_pos, remark1, remark2, &
                                                title
    INTEGER                                  :: iskip, nat, outformat, &
                                                traj_unit
    LOGICAL                                  :: failure, init_traj
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    NULLIFY(logger, cell, subsys, particles, particle_set)
    logger => cp_error_get_logger(error)
    my_pos = "APPEND"
    my_act = "WRITE"

    CALL force_env_get(force_env, cell=cell, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys=subsys, particles=particles,error=error)
    particle_set => particles%els
    nat=particles%n_els

    CALL get_output_format(root_section,"MOTION%PRINT%TRAJECTORY",my_form,my_ext,error)
    traj_unit = cp_print_key_unit_nr(logger,root_section,"MOTION%PRINT%TRAJECTORY",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="pos", is_new_file=init_traj,&
            error=error)
    IF (traj_unit > 0) THEN
       CALL section_vals_val_get(root_section,"MOTION%PRINT%TRAJECTORY%FORMAT",i_val=outformat,error=error)
       title = ""
       IF     (outformat==dump_dcd) THEN
          IF ( init_traj ) THEN
             !Lets write the header for the coordinate dcd
             iskip = section_get_ival(root_section,"MOTION%PRINT%TRAJECTORY%EACH",error=error)
             WRITE(traj_unit) "CORD",0,it,iskip,0,0,0,0,0,0,REAL(it,KIND=sp),1,0,0,0,0,0,0,0,0,24
             remark1= "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
             remark2= "REMARK Support new DCD format with cell information"
             WRITE(traj_unit) 2,remark1,remark2
             WRITE(traj_unit) nat
             CALL m_flush(traj_unit)
          ENDIF
       ELSE IF(outformat==dump_xmol) THEN
          !         *** This information can be digested by Molden ***
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i = ",it," , E =",etot
       ELSE IF(outformat==dump_atomic) THEN
          ! Do Nothing..
       ELSE
          CPPostcondition(.FALSE.,cp_fatal_level,routineP,error,failure)
       END IF
       CALL write_particle_coordinates(particle_set,traj_unit,outformat,"POS",TRIM(title),cell)
    END IF
    CALL cp_print_key_finished_output(traj_unit,logger,root_section,"MOTION%PRINT%TRAJECTORY",error=error)
  END SUBROUTINE write_geo_traj

END MODULE gopt_f_methods
