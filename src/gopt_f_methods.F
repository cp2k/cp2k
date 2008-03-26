!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains a functional that calculates the energy and its derivatives
!>      for the geometry optimizer
!> \par History
!>      none
! *****************************************************************************
MODULE gopt_f_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE dimer_types,                     ONLY: dimer_env_type
  USE dimer_utils,                     ONLY: update_dimer_vec
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_cell_method_id,&
                                             default_minimization_method_id,&
                                             default_ts_method_id
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE motion_utils,                    ONLY: write_simulation_cell,&
                                             write_stress_tensor,&
                                             write_trajectory
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_structure_data
  USE qmmm_util,                       ONLY: apply_qmmm_translate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_evaluate
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_f_methods'

  PUBLIC :: gopt_f_create_x0, write_cycle_infos, check_converg, write_final_info,&
       print_geo_opt_header, check_rot_conv, write_rot_cycle_infos, print_geo_opt_nc,&
       geo_opt_io, gopt_f_io_init, gopt_f_io, gopt_f_io_finalize, gopt_f_ii

CONTAINS

! *****************************************************************************
!> \brief returns the value of the parameters for the actual configuration
!> \param gopt_env the geometry optimization environment you want the info about
!>      x0: the parameter vector (is allocated by this routine)
!> \param error variable to control error logging, stopping,... 
!>             see module cp_error_handling 
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE gopt_f_create_x0(gopt_env, x0, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_create_x0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idg, idir, ip, j, ndg, stat
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    NULLIFY(particles,subsys,cell)

    SELECT CASE(gopt_env%type_id)
    CASE(default_minimization_method_id,default_ts_method_id)
       CALL force_env_get(gopt_env%force_env, subsys=subsys,error=error)
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
    CASE(default_cell_method_id)
       CALL force_env_get(gopt_env%force_env, cell=cell, error=error)
       ALLOCATE(x0(6),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       idg = 0
       DO i = 1, 3
          DO j = 1, i
             idg = idg + 1
             x0(idg) = cell%hmat(j,i)
          END DO
       END DO
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE gopt_f_create_x0

! *****************************************************************************
!> \brief Prints iteration step of the optimization procedure on screen
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE gopt_f_ii(its, output_unit)

    INTEGER, INTENT(IN)                      :: its, output_unit

    IF ( output_unit > 0) THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,26('-'))")
       WRITE(UNIT=output_unit,FMT="(T2,A,I6)")"OPTIMIZATION STEP: ",its
       WRITE(UNIT=output_unit,FMT="(T2,26('-'))")
       CALL m_flush(output_unit)
    END IF
  END SUBROUTINE gopt_f_ii

! *****************************************************************************
!> \brief Handles the Output during an optimization run 
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE gopt_f_io_init(gopt_env, output_unit, opt_energy, wildcard, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(KIND=dp)                            :: opt_energy
    CHARACTER(LEN=5)                         :: wildcard
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_io_init', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: pres_int

    SELECT CASE(gopt_env%type_id)
    CASE(default_ts_method_id,default_minimization_method_id)
       ! Geometry Optimization (Minimization and Transition State Search)
       IF (.NOT.gopt_env%dimer_rotation) THEN
          CALL write_cycle_infos(output_unit,it=0,etot=opt_energy,wildcard=wildcard)
       ELSE
          CALL write_rot_cycle_infos(output_unit,it=0,etot=opt_energy,dimer_env=gopt_env%dimer_env,&
               wildcard=wildcard)
       END IF
    CASE(default_cell_method_id)
       ! Cell Optimization
       pres_int = gopt_env%cell_env%pres_int
       CALL write_cycle_infos(output_unit,it=0,etot=opt_energy,pres_int=pres_int,wildcard=wildcard)
    END SELECT

  END SUBROUTINE gopt_f_io_init

! *****************************************************************************
!> \brief Handles the Output during an optimization run 
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE gopt_f_io(gopt_env, force_env, root_section, its, opt_energy,&
       output_unit, eold, emin, wildcard, gopt_param, ndf, dx, xi, conv, ediff, pred, rat,&
       step, rad, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: its
    REAL(KIND=dp), INTENT(IN)                :: opt_energy
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(KIND=dp)                            :: eold, emin
    CHARACTER(LEN=5)                         :: wildcard
    TYPE(gopt_param_type), POINTER           :: gopt_param
    INTEGER, INTENT(IN), OPTIONAL            :: ndf
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: dx
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: xi
    LOGICAL, OPTIONAL                        :: conv
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: ediff, pred, rat, step, rad
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_io', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pres_diff, pres_int, pres_tol

    failure = .FALSE.
    SELECT CASE(gopt_env%type_id)
    CASE(default_ts_method_id,default_minimization_method_id)
       ! Geometry Optimization (Minimization and Transition State Search)
       IF (.NOT.gopt_env%dimer_rotation) THEN
          CALL geo_opt_io(force_env=force_env, root_section=root_section, &
               motion_section=gopt_env%motion_section, its=its, opt_energy=opt_energy, error=error)
          CALL write_cycle_infos(output_unit,its,etot=opt_energy,ediff=opt_energy-eold,&
               pred=pred,rat=rat,step=step,rad=rad,emin=emin,wildcard=wildcard)
          ! Possibly check convergence
          IF (PRESENT(conv)) THEN 
             CPPostcondition(PRESENT(ndf),cp_failure_level,routineP,error,failure)
             CPPostcondition(PRESENT(dx),cp_failure_level,routineP,error,failure)
             CPPostcondition(PRESENT(xi),cp_failure_level,routineP,error,failure)
             CALL check_converg(ndf,dx,xi,output_unit,conv,gopt_param)
          END IF
       ELSE
          CALL update_dimer_vec(gopt_env%dimer_env, gopt_env%motion_section, error)
          CALL write_restart(force_env=force_env,root_section=root_section,&
               error=error)
          CALL write_rot_cycle_infos(output_unit,its,opt_energy,opt_energy-eold,emin,gopt_env%dimer_env,&
               wildcard)
          ! Possibly check convergence
          IF (PRESENT(conv)) THEN 
             CPPostcondition(ASSOCIATED(gopt_env%dimer_env),cp_failure_level,routineP,error,failure)
             CALL check_rot_conv(gopt_env%dimer_env, output_unit, conv, error)
          END IF
       END IF
    CASE(default_cell_method_id)
       ! Cell Optimization
       pres_diff = gopt_env%cell_env%pres_int - gopt_env%cell_env%pres_ext
       pres_int  = gopt_env%cell_env%pres_int
       pres_tol  = gopt_env%cell_env%pres_tol
       CALL geo_opt_io(force_env=force_env, root_section=root_section, &
            motion_section=gopt_env%motion_section, its=its, opt_energy=opt_energy, error=error)
       CALL write_cycle_infos(output_unit,its,etot=opt_energy,ediff=opt_energy-eold,&
            pred=pred,rat=rat,step=step,rad=rad,emin=emin,pres_int=pres_int,wildcard=wildcard)
       ! Possibly check convergence
       IF (PRESENT(conv)) THEN 
          CPPostcondition(PRESENT(ndf),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(dx),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(xi),cp_failure_level,routineP,error,failure)
          CALL check_converg(ndf,dx,xi,output_unit,conv,gopt_param,pres_diff,pres_tol)
       END IF
    END SELECT
  END SUBROUTINE gopt_f_io

! *****************************************************************************
!> \brief Handles the Output at the end of an optimization run 
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE gopt_f_io_finalize(gopt_env, force_env, x0, conv, its, root_section,&
       para_env, master, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    LOGICAL                                  :: conv
    INTEGER                                  :: its
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: master
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_io_finalize', &
      routineP = moduleN//':'//routineN

    IF (gopt_env%eval_opt_geo) THEN
       IF (.NOT.gopt_env%dimer_rotation) THEN
          CALL write_final_info(conv, its, gopt_env, x0, master, para_env, force_env,&
               gopt_env%motion_section, root_section, error)
       ELSE
          CALL update_dimer_vec(gopt_env%dimer_env, gopt_env%motion_section, error)
          CALL write_restart(force_env=force_env,root_section=root_section,&
               error=error)     
       END IF
    END IF

  END SUBROUTINE gopt_f_io_finalize

! *****************************************************************************
  SUBROUTINE write_cycle_infos(output_unit,it,etot,ediff,pred,rat,step,rad,emin,&
       pres_int,wildcard)

    INTEGER, INTENT(IN)                      :: output_unit, it
    REAL(KIND=dp), INTENT(IN)                :: etot
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: ediff, pred, rat, step, rad, &
                                                emin, pres_int
    CHARACTER(LEN=5), OPTIONAL               :: wildcard

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_cycle_infos', &
      routineP = moduleN//':'//routineN

    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(unit)
    CALL get_cp2k_units(unit)
    IF (output_unit>0)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,8('-'),A,I5,1X,12('-'))")&
            "  Informations at step = ",it
       WRITE(UNIT=output_unit,FMT="(T2,A,T47,A)")&
            " Optimization Method        = ",wildcard
       WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
            " Total Energy               = ",etot
       IF (PRESENT(pres_int)) THEN
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
               " Internal Pressure [bar]    = ",pres_int*unit%bar
       END IF
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
    END IF
  END SUBROUTINE write_cycle_infos

! *****************************************************************************
!> \author Luca Bellucci and Teodoro Laino - created [tlaino]
!> \data  01.2008
! *****************************************************************************
  SUBROUTINE write_rot_cycle_infos(output_unit,it,etot,ediff,emin,dimer_env,wildcard)

    INTEGER, INTENT(IN)                      :: output_unit, it
    REAL(KIND=dp), INTENT(IN)                :: etot
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: ediff, emin
    TYPE(dimer_env_type), POINTER            :: dimer_env
    CHARACTER(LEN=5), OPTIONAL               :: wildcard

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_rot_cycle_infos', &
      routineP = moduleN//':'//routineN

    IF (output_unit>0)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,4('-'),A,I5,1X,5('-'))")&
            "  Informations at rotational step = ",it
       WRITE(UNIT=output_unit,FMT="(T2,A,T47,A)")&
            " Optimization Method        = ",wildcard
       WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
            " Local Curvature            = ",dimer_env%rot%curvature
       WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
            " Total Rotational Force     = ",etot
       IF (PRESENT(ediff)) THEN
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
               " Real Force change          = ",ediff
       END IF
       IF (PRESENT(emin)) THEN
          IF( etot < emin )THEN
             WRITE(UNIT=output_unit,FMT="(T2,2A)")&
                  " Decrease in rotational force =   ",&
                  "             YES"
          ELSE
             WRITE(UNIT=output_unit,FMT="(T2,2A)")&
                  " Decrease in rotational force =   ",&
                  "              NO"
          END IF
       END IF
       IF (it==0) WRITE(UNIT=output_unit,FMT="(T2,51('-'))")
    END IF
  END SUBROUTINE write_rot_cycle_infos

! *****************************************************************************
  SUBROUTINE check_converg(ndf,dr,g,output_unit,conv,gopt_param,pres_diff,pres_tol)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(IN)                :: dr(ndf), g(ndf)
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(OUT)                     :: conv
    TYPE(gopt_param_type), POINTER           :: gopt_param
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: pres_diff, pres_tol

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_converg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: indf
    LOGICAL                                  :: conv_dx, conv_g, conv_p, &
                                                conv_rdx, conv_rg
    REAL(KIND=dp)                            :: dumm, dxcon, gcon, maxdum(4), &
                                                rmsgcon, rmsxcon
    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(unit)
    CALL get_cp2k_units(unit)
    dxcon = gopt_param%max_dr
    gcon = gopt_param%max_force
    rmsgcon = gopt_param%rms_force
    rmsxcon = gopt_param%rms_dr

    conv     = .FALSE.
    conv_dx  = .TRUE.
    conv_rdx = .TRUE.
    conv_g   = .TRUE.
    conv_rg  = .TRUE.
    conv_p   = .TRUE.

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

    IF (PRESENT(pres_diff).AND.PRESENT(pres_tol)) conv_p = ABS(pres_diff)<ABS(pres_tol)

    IF(output_unit>0)THEN
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
       IF (PRESENT(pres_diff).AND.PRESENT(pres_tol)) THEN
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
               " Pressure Deviation [bar]   = ",pres_diff*unit%bar
          WRITE(UNIT=output_unit,FMT="(T2,A,F20.10)")&
               " Pressure Tolerance [bar]   = ",pres_tol*unit%bar
          IF(conv_p)THEN
             WRITE(UNIT=output_unit,FMT="(T2,2A)")&
                  " Conv. for  PRESSURE        =     ",&
                  "             YES"
             
          ELSE
             WRITE(UNIT=output_unit,FMT="(T2,2A)")&
                  " Conv. for  PRESSURE        =     ",&
                  "              NO"
          END IF
       END IF
       WRITE(UNIT=output_unit,FMT="(T2,51('-'))")
    END IF

    IF (conv_dx .AND. conv_rdx .AND. conv_g .AND. conv_rg .AND. conv_p) conv = .TRUE.

    IF ((conv).AND.(output_unit>0)) THEN
       WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       WRITE (UNIT=output_unit,FMT="(T2,A,T25,A,T78,A)")&
            "***","GEOMETRY OPTIMIZATION COMPLETED","***"
       WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
    END IF

  END SUBROUTINE check_converg

! *****************************************************************************
!> \author Luca Bellucci and Teodoro Laino - created [tlaino]
!> \data  01.2008
! *****************************************************************************
  SUBROUTINE check_rot_conv(dimer_env, output_unit, conv, error)

    TYPE(dimer_env_type), POINTER            :: dimer_env
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(OUT)                     :: conv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_rot_conv', &
      routineP = moduleN//':'//routineN

    conv = (ABS(dimer_env%rot%angle2)< dimer_env%rot%angle_tol) 
    IF(output_unit>0)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
            " Convergence check :"
       WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
            " Predicted angle step size  =     ",dimer_env%rot%angle1
       WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
            " Effective angle step size  =     ",dimer_env%rot%angle2
       WRITE(UNIT=output_unit,FMT="(T2,A,F16.10)")&
            " Conv. limit for angle step size =",dimer_env%rot%angle_tol
       IF(conv)THEN
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
               " Convergence in angle step size   =",&
               "            YES"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,2A)")&
               " Convergence in angle step size   =",&
               "             NO"
       END IF
       WRITE(UNIT=output_unit,FMT="(T2,51('-'))")
    END IF
    IF ((conv).AND.(output_unit>0)) THEN
       WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       WRITE (UNIT=output_unit,FMT="(T2,A,T25,A,T78,A)")&
            "***","ROTATION OPTIMIZATION COMPLETED","***"
       WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
    END IF

  END SUBROUTINE check_rot_conv

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich 
!> \data  11.2007
! *****************************************************************************
  RECURSIVE SUBROUTINE write_final_info(conv, it, gopt_env, x0, master, para_env, force_env,&
    motion_section, root_section, error)
    LOGICAL, INTENT(IN)                      :: conv
    INTEGER, INTENT(INOUT)                   :: it
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    INTEGER, INTENT(IN)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: motion_section, root_section
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
                         error=error)

      ! Write out more data for the minimum
      IF (logger%para_env%mepos==logger%para_env%source) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger), fmt="(/,T20,' Reevaluating energy at the minimum ')")
      END IF
      CALL cp_eval_at(gopt_env,x0,f=etot,master=master,para_env=para_env,error=error)
      CALL write_geo_traj(force_env,  root_section, it, etot, error)
    END IF

  END SUBROUTINE write_final_info

! *****************************************************************************
!> \brief  Specific driver for dumping trajectory during a GEO_OPT
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \data   11.2007
! *****************************************************************************
  SUBROUTINE write_geo_traj(force_env,  root_section, it, etot, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: it
    REAL(KIND=dp), INTENT(IN)                :: etot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_geo_traj', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL write_trajectory(force_env,  root_section, it, 0.0_dp, etot, error=error)

  END SUBROUTINE write_geo_traj

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \data  01.2008
! *****************************************************************************
  SUBROUTINE print_geo_opt_header(gopt_env, output_unit, label)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    INTEGER, INTENT(IN)                      :: output_unit
    CHARACTER(LEN=*), INTENT(IN)             :: label

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_geo_opt_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_format, my_label
    INTEGER                                  :: ix

    IF(output_unit>0)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       IF (gopt_env%dimer_rotation) THEN
          my_label = "OPTIMIZING DIMER ROTATION"
       ELSE
          my_label = "STARTING "//gopt_env%tag(1:8)//" OPTIMIZATION"
       END IF

       ix = (80 - 7 - LEN_TRIM(my_label))/2
       ix = ix + 5
       my_format="(T2,A,T"//cp_to_string(ix)//",A,T78,A)"
       WRITE(UNIT=output_unit,FMT=TRIM(my_format))"***",TRIM(my_label),"***"

       ix = (80 - 7 - LEN_TRIM(label))/2
       ix = ix + 5
       my_format="(T2,A,T"//cp_to_string(ix)//",A,T78,A)"
       WRITE(UNIT=output_unit,FMT=TRIM(my_format))"***",TRIM(label),"***"

       WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       CALL m_flush(output_unit)
    END IF
  END SUBROUTINE print_geo_opt_header

! *****************************************************************************
!> \author Teodoro Laino [tlaino] - University of Zurich 
!> \data  01.2008
! *****************************************************************************
  SUBROUTINE print_geo_opt_nc(gopt_env,output_unit)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    INTEGER, INTENT(IN)                      :: output_unit

    IF(output_unit>0)THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
            "*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***"
       IF (.NOT.gopt_env%dimer_rotation) THEN
          WRITE(UNIT=output_unit,FMT="(T2,A)")&
               "***        EXITING GEOMETRY OPTIMIZATION         ***"
       ELSE
          WRITE(UNIT=output_unit,FMT="(T2,A)")&
               "***        EXITING ROTATION OPTIMIZATION         ***"
       END IF
       CALL m_flush(output_unit)
    END IF
  END SUBROUTINE print_geo_opt_nc

! *****************************************************************************
!> \brief   Prints information during GEO_OPT common to all optimizers
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    02.2008
!> \version 1.0
! *****************************************************************************
  SUBROUTINE geo_opt_io(force_env, root_section, motion_section, its, opt_energy, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section, motion_section
    INTEGER, INTENT(IN)                      :: its
    REAL(KIND=dp), INTENT(IN)                :: opt_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geo_opt_io', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(unit, para_env, atomic_kind_set, subsys, particle_set, local_particles,&
            atomic_kinds, particles)
    CALL get_cp2k_units(unit)

    ! Write Restart File
    CALL write_restart(force_env=force_env,root_section=root_section,&
         error=error)

    ! Write Trajectory
    CALL write_geo_traj(force_env,  root_section, its, opt_energy, error)

    ! Write the stress Tensor
    CALL force_env_get(force_env, cell=cell, virial=virial, para_env=para_env, &
         subsys=subsys, error=error)
    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles, error=error)
    atomic_kind_set => atomic_kinds%els
    particle_set => particles%els
    CALL virial_evaluate(atomic_kind_set, particle_set,  local_particles, &
         virial, para_env%group, error=error)
    CALL write_stress_tensor(virial, cell, motion_section, unit, its, 0.0_dp, error=error)
    
    ! Write the cell
    CALL write_simulation_cell(cell, motion_section, its, 0.0_dp, error=error)

  END SUBROUTINE geo_opt_io

END MODULE gopt_f_methods
