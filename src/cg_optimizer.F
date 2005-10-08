!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/cg_optimizer [1.0] *
!!
!!   NAME
!!     cg_optimizer
!!
!!   FUNCTION
!!     Routines for Geometry optimization using  Conjugate Gradients
!!
!!   AUTHOR
!!     Teodoro Laino [teo]
!!     10.2005
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************

MODULE cg_optimizer
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_param_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_particle_coordinates,&
                                             write_structure_data
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE bfgs_optimizer,                  ONLY: write_cycle_infos,&
                                             check_converg
  USE input_constants,                 ONLY: dump_xmol
  USE cp_files,                        ONLY: close_file,open_file
  USE memory_utilities,                ONLY: reallocate

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: geoopt_cg

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cg_optimizer'

CONTAINS

!!****f* cg_optimizer/geoopt_cg [1.0] *
!!
!!   NAME
!!     geoopt_cg
!!
!!   FUNCTION
!!     Driver for conjugate gradient optimization technique
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** ********************************************************************** 
  SUBROUTINE geoopt_cg(force_env,gopt_param,globenv,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'geoopt_cg', &
      routineP = moduleN//':'//routineN
    INTEGER :: atom, handle, iatom, ikind, indf, j, &
      natom, ndf, nfixed, nfree, nkind, output_unit, stat, dataunit
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: conv, ionode
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL write_checkpoint_information("entering "//routineN,globenv)
    CALL timeset(routineN,"I","",handle)   
    IF (.NOT.failure) THEN
       ionode      = globenv%ionode
       output_unit = globenv%scr

       conv        = .FALSE.
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
       nfixed = COUNT(particle_set(:)%is_fixed)
       nfree = natom - nfixed
       ndf = 3*nfree
       
       ALLOCATE(x(ndf), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
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
                END DO
             END IF
          END DO
       END DO
       
       IF(ionode)THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T28,A,T78,A)")&
               "***","STARTING GEOMETRY OPTIMIZATION","***"
          WRITE(UNIT=output_unit,FMT="(T2,A,T28,A,T78,A)")&
               "***","    CONJUGATE GRADIENTS","***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(/,(T2,A,I10))")&
               "Number of free atoms         = ",nfree,&
               "Number of degrees of freedom = ",ndf
          CALL m_flush(output_unit)
          
          CALL open_file(file_name="bfgs.data",&
                         file_status="UNKNOWN",&
                         file_form="FORMATTED",&
                         file_action="WRITE",&
                         unit_number=dataunit)
       END IF

       CALL write_qs_particle_coordinates(particle_set,cell,globenv) 
       CALL write_particle_distances(particle_set,cell,globenv)
       CALL write_structure_data(particle_set,cell,globenv)
       
       CALL  cp_cg_main(force_env, x, gopt_param, ionode, output_unit, dataunit, &
            globenv, error)

       CALL write_particle_distances(particle_set,cell,globenv)
       CALL write_structure_data(particle_set,cell,globenv)
   
       DEALLOCATE(x, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END SUBROUTINE geoopt_cg

!!****f* cg_optimizer/get_energy_and_gradients [1.0] *
!!
!!   NAME
!!     get_energy_and_gradients
!!
!!   FUNCTION
!!     Main driver to evaluate energy and derivatives
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  SUBROUTINE get_energy_and_gradients(force_env,xvec,energy,grad)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: energy
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: grad
    ! Local Variables
    INTEGER                                  :: atom, iatom, ikind, indf, j, &
                                                natom, nkind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: get_energy, get_gradients
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    get_gradients = PRESENT(grad)
    get_energy    = PRESENT(energy)
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
                particle_set(atom)%r(j) = xvec(indf)
             END DO
          END IF
       END DO
    END DO

    CALL force_env_calc_energy_force(force_env,calc_force=get_gradients)

    IF (get_gradients)  THEN
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
                   grad(indf) = -particle_set(atom)%f(j)
                END DO
             END IF
          END DO
       END DO
    END IF
    IF (get_energy) CALL force_env_get(force_env=force_env,potential_energy=energy)    
  END SUBROUTINE get_energy_and_gradients

!!****f* cg_optimizer/cp_cg_main [1.0] *
!!
!!   NAME
!!     get_energy_and_gradients
!!
!!   FUNCTION
!!     This really performs the conjugate gradients optimization
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  SUBROUTINE cp_cg_main(force_env, xvec, gopt_param, ionode, output_unit, dataunit,&
       globenv, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec
    INTEGER, INTENT(IN)                      :: output_unit, dataunit
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cg_main', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: its, stat, maxiter, brent_max_iter
    REAL(KIND=dp)                            :: dgg, gam, gg, emin, eold, opt_energy
    REAL(KIND=dp)                            :: brent_tol, brack_limit, step, res_lim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g, h, xi
    LOGICAL                                  :: conv
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: xold
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    CHARACTER(LEN=80)                        :: title
    CHARACTER(LEN=5)                         :: wildcard
    LOGICAL                                  :: failure,  Fletcher_Reeves_ON, use_fit
    INTEGER                                  :: max_steep_steps

    failure = .FALSE.
    NULLIFY(g, h, xi)
    IF (.NOT.failure) THEN
       maxiter            = gopt_param%max_iter
       max_steep_steps    = gopt_param%max_steep_steps
       Fletcher_Reeves_ON = gopt_param%Fletcher_Reeves_ON
       brent_tol          = gopt_param%brent_tol
       brent_max_iter     = gopt_param%brent_max_iter
       brack_limit        = gopt_param%brack_limit
       use_fit            = gopt_param%use_fit_line_search
       res_lim            = gopt_param%restart_limit
       ALLOCATE(g (SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(h (SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xi(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xold(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL get_energy_and_gradients(force_env, xvec, opt_energy, xi)
       g  = -xi
       h  =   g
       xi =   h
       emin    = HUGE(0.0_dp)
       ! Main Loop
       wildcard = "   SD"
       DO its = 1, maxiter
          xold = xvec
          eold = opt_energy
          emin = MIN(emin, opt_energy)
          step = 1.0_dp
          IF (use_fit) THEN
             CALL linminE(force_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, brack_limit, &
                          step, output_unit, ionode, gopt_param, error)             
          ELSE
             CALL linmin(force_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, brack_limit, &
                         step, output_unit, ionode, error)
          END IF
          !
          ! Convergence check
          !       
          CALL write_cycle_infos(output_unit,its,opt_energy,opt_energy-eold,emin=emin,wildcard=wildcard)
          CALL check_converg(SIZE(xvec),xvec-xold,xi,ionode,output_unit,conv,gopt_param)
          wildcard = "   CG"
          IF (ionode) THEN
             CALL force_env_get(force_env,&
                                cell=cell,&
                                subsys=subsys)           
             CALL cp_subsys_get(subsys=subsys(1)%subsys,&
                                particles=particles)
             particle_set => particles%els              

             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)") " i =",its,", E =",opt_energy
             CALL write_particle_coordinates(particle_set,cell,dataunit,&
                                             dump_xmol,"POS",TRIM(title))
             CALL m_flush(dataunit)
             CALL write_qs_particle_coordinates(particle_set,cell,globenv)
          END IF
          IF (conv.OR.(its==maxiter))  EXIT
          CALL get_energy_and_gradients(force_env, xvec, opt_energy, xi)          
          gg  = DOT_PRODUCT(g,g)
          IF (Fletcher_Reeves_ON) THEN 
             dgg = DOT_PRODUCT(xi,xi)
          ELSE
             dgg = DOT_PRODUCT((xi+g),xi)
          END IF
          gam = dgg/gg
          g   = -xi
          h   = g + gam * h
          xi  = h
          !
          ! Reset Condition or Steepest Descent Requested
          !
          IF ( ABS(DOT_PRODUCT(g,h))/SQRT((DOT_PRODUCT(g,g)*DOT_PRODUCT(h,h))) < res_lim &
               .OR. its+1 <= max_steep_steps) THEN
             g  = -xi
             h  =   g
             xi =   h
             wildcard = "   SD"
          END IF
       END DO
       IF(its == maxiter .AND. (.NOT.conv))THEN
          IF(ionode)THEN
             WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                  "*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***"
             WRITE(UNIT=output_unit,FMT="(T2,A)")&
                  "***        EXITING GEOMETRY OPTIMIZATION         ***" 
             CALL m_flush(output_unit)
          END IF
       END IF
       DEALLOCATE(xold, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(g, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(h, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(xi,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE cp_cg_main

!!****f* cg_optimizer/LinMinE [1.0] *
!!
!!   NAME
!!     LinMinE
!!
!!   FUNCTION
!!     Line Minimization routine Enhanced
!!
!!   NOTES
!!     Given as input the vector XVEC and XI, finds the scalar
!!     xmin that minimizes the energy XVEC+xmin*XI. Replaces XMIN
!!     with the optimal value. Enhanced Version
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  SUBROUTINE LinMinE(force_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, &
       brack_limit, step, output_unit, ionode, gopt_param, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, brent_tol, brack_limit, step
    INTEGER                                  :: brent_max_iter, output_unit
    LOGICAL                                  :: ionode
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'LinMinE', &
      routineP = moduleN//':'//routineN    
    INTEGER                                  :: stat, odim, loc_iter
    REAL(KIND=dp)                            :: ax, bx, xmin, xx, fprev, rms_force, rms_dr, scale
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: hist
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(pcom,xicom,hist)
    rms_dr    = gopt_param%rms_dr
    rms_force = gopt_param%rms_force
    IF (.NOT.failure) THEN
       ALLOCATE(pcom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xicom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       pcom  = xvec
       xi = xi / SQRT(DOT_PRODUCT(xi,xi))
       step = 0.05_dp
       xicom = xi
       ax = 0.0_dp
       xx = step
       CALL cg_mnbrak(force_env, ax, xx, bx, pcom, xicom, brack_limit, output_unit, ionode, &
            histpoint=hist,error=error)
       !
       fprev      = 0.0_dp
       opt_energy = MINVAL(hist(:,2))
       odim =  SIZE(hist,1)
       scale = 0.25_dp
       loc_iter = 0
       DO WHILE ( ABS(hist(odim,3)) > rms_force*scale .OR. ABS(hist(odim,1)-hist(odim-1,1)) > scale*rms_dr)
          loc_iter = loc_iter + 1
          fprev = opt_energy
          xmin  = FindMin(hist(:,1),hist(:,2),hist(:,3),error)
          CALL reallocate(hist,1,odim+1,1,3)
          hist(odim+1,1) = xmin
          hist(odim+1,3) = cg_deval1d(force_env,xmin,pcom,xicom,opt_energy,error)
          hist(odim+1,2) = opt_energy
          odim = SIZE(hist,1)
       END DO
       !
       xi = xmin*xi
       xvec = xvec + xi
       DEALLOCATE(pcom,  stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(xicom, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(hist, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ionode) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
               "***","FIT LS  - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       END IF       
    END IF
  END SUBROUTINE LinMinE

!!****f* cg_optimizer/LinMin [1.0] *
!!
!!   NAME
!!     LinMin
!!
!!   FUNCTION
!!     Line Minimization routine
!!
!!   NOTES
!!     Given as input the vector XVEC and XI, finds the scalar
!!     xmin that minimizes the energy XVEC+xmin*XI. Replaces XMIN
!!     with the optimal value
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  SUBROUTINE LinMin(force_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, &
       brack_limit, step, output_unit, ionode, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, brent_tol, brack_limit, step
    INTEGER                                  :: brent_max_iter, output_unit
    LOGICAL                                  :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'LinMin', &
      routineP = moduleN//':'//routineN    
    INTEGER                                  :: stat
    REAL(KIND=dp)                            :: ax, bx, xmin, xx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(pcom,xicom)
    IF (.NOT.failure) THEN
       ALLOCATE(pcom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xicom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       pcom  = xvec
       xicom = xi
       ax = 0.0_dp
       xx = step       
       CALL cg_mnbrak(force_env, ax, xx, bx, pcom, xicom, brack_limit, output_unit, ionode,&
            error=error)
       opt_energy = cg_dbrent(force_env, ax, xx, bx, brent_tol, brent_max_iter, &
                              xmin, pcom, xicom, output_unit, ionode, error)
       xi = xmin*xi
       xvec = xvec + xi
       DEALLOCATE(pcom,  stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(xicom, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE LinMin

!!****f* cg_optimizer/cg_mnbrak [1.0] *
!!
!!   NAME
!!     cg_mnbrak
!!
!!   FUNCTION
!!     Routine for intially bracketing a minimum based on the golden search
!!     minimum
!!
!!   NOTES
!!     Given two distinct initial points ax and bx this routine searches 
!!     in the downhill direction and returns new points ax, bx, cx that
!!     bracket the minimum of the function 
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  SUBROUTINE cg_mnbrak(force_env, ax, bx, cx, pcom, xicom, brack_limit, output_unit, ionode,&
       histpoint, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp)                            :: ax, bx, cx, brack_limit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp), DIMENSION(:,:), POINTER,&
         OPTIONAL                            :: histpoint
    INTEGER                                  :: output_unit
    LOGICAL                                  :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cg_mnbrak', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                            :: dum, fa, fb, fc, fu, q, r, u, &
                                                ulim, gold
    LOGICAL                                  :: failure, hist
    INTEGER                                  :: loc_iter, odim, stat

    failure = .FALSE.
    hist    = PRESENT(histpoint)
    IF (hist) THEN
       CPPrecondition(.NOT.ASSOCIATED(histpoint),cp_failure_level,routineP,error,failure)
       ALLOCATE(histpoint(3,3), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    gold = (1.0_dp+SQRT(5.0_dp))/2.0_dp 
    IF (.NOT.failure) THEN
       IF (hist) THEN
          histpoint(1,1) = ax
          histpoint(1,3) = cg_deval1d(force_env,ax,pcom,xicom,fa,error)
          histpoint(1,2) = fa
          histpoint(2,1) = bx
          histpoint(2,3) = cg_deval1d(force_env,bx,pcom,xicom,fb,error)
          histpoint(2,2) = fb
       ELSE
          fa=cg_eval1d(force_env,ax,pcom,xicom) 
          fb=cg_eval1d(force_env,bx,pcom,xicom) 
       END IF
       IF(fb.GT.fa)THEN
          dum=ax
          ax=bx
          bx=dum
          dum=fb
          fb=fa
          fa=dum
       ENDIF
       cx=bx+gold*(bx-ax)
       IF (hist) THEN
          histpoint(3,1) = cx
          histpoint(3,3) = cg_deval1d(force_env,cx,pcom,xicom,fc,error)
          histpoint(3,2) = fc          
       ELSE
          fc=cg_eval1d(force_env,cx,pcom,xicom)
       END IF
       loc_iter = 3
       DO WHILE (fb.GE.fc)
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*SIGN(MAX(ABS(q-r),TINY(0.0_dp)),q-r))
          ulim=bx+brack_limit*(cx-bx)
          IF ((bx-u)*(u-cx).GT.0.0_dp) THEN
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(force_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(force_env,u,pcom,xicom)
             END IF
             loc_iter = loc_iter + 1
             IF(fu.LT.fc)THEN
                ax=bx
                fa=fb
                bx=u
                fb=fu
                EXIT
             ELSE IF(fu.GT.fb)THEN
                cx=u
                fc=fu
                EXIT
             ENDIF
             u=cx+gold*(cx-bx)
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(force_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(force_env,u,pcom,xicom)
             END IF
             loc_iter = loc_iter + 1
          ELSE IF((cx-u)*(u-ulim).GT.0.)THEN
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(force_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(force_env,u,pcom,xicom)
             END IF
             loc_iter = loc_iter + 1
             IF(fu.LT.fc)THEN
                bx=cx
                cx=u
                u=cx+gold*(cx-bx)
                fb=fc
                fc=fu
                IF (hist) THEN
                   odim = SIZE(histpoint,1)
                   CALL reallocate(histpoint,1,odim+1,1,3)
                   histpoint(odim+1,1) = u
                   histpoint(odim+1,3) = cg_deval1d(force_env,u,pcom,xicom,fu,error)
                   histpoint(odim+1,2) = fu
                ELSE
                   fu=cg_eval1d(force_env,u,pcom,xicom)
                END IF
                loc_iter = loc_iter + 1
             ENDIF
          ELSE IF((u-ulim)*(ulim-cx).GE.0.)THEN
             u=ulim
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(force_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(force_env,u,pcom,xicom)
             END IF
             loc_iter = loc_iter + 1
          ELSE
             u=cx+gold*(cx-bx)
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(force_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(force_env,u,pcom,xicom)
             END IF
             loc_iter = loc_iter + 1
          ENDIF
          ax=bx
          bx=cx
          cx=u
          fa=fb
          fb=fc
          fc=fu
       END DO
       IF (ionode) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
               "***","MNBRACK - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       END IF
    END IF
  END SUBROUTINE cg_mnbrak

!!****f* cg_optimizer/cg_dbrent [1.0] *
!!
!!   NAME
!!     cg_dbrent
!!
!!   FUNCTION
!!     Routine implementing the Brent Method
!!     Brent,R.P. Algorithm for Minimization without Derivatives, Chapt.5
!!     1973
!!     Extension in the use of derivatives
!!
!!   NOTES
!!     Given a bracketing  triplet of abscissas ax, bx, cx (such that bx 
!!     is between ax and cx and energy of bx is less than energy of ax and cx),
!!     this routine isolates the minimum to a precision of about tol using 
!!     Brent method. This routine implements the extension of the Brent Method
!!     using derivatives
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  FUNCTION cg_dbrent(force_env,ax,bx,cx,tol,itmax,xmin,pcom,xicom,output_unit, ionode, error) RESULT(dbrent)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp)                            :: ax, bx, cx, tol, xmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: dbrent
    INTEGER                                  :: itmax
    INTEGER                                  :: output_unit
    LOGICAL                                  :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cg_dbrent', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: zeps = 1.0E-8_dp
    INTEGER                                  :: iter, loc_iter
    LOGICAL                                  :: ok1, ok2, skip0, skip1
    REAL(KIND=dp)                            :: a, b, d, d1, d2, du, dv, dw, &
                                                dx, e, fu, fv, fw, fx, olde, &
                                                tol1, tol2, u, u1, u2, v, w, &
                                                x, xm
    LOGICAL                                  :: failure
    
    failure = .FALSE.
    IF (.NOT.failure) THEN
       a  = MIN(ax,cx)
       b  = MAX(ax,cx)
       v  = bx ; w  = v; x  = v
       e  = 0.0_dp
       dx = cg_deval1d(force_env,x,pcom,xicom,fx,error)
       fv = fx
       fw = fx
       dv = dx
       dw = dx
       loc_iter = 1
       DO  iter = 1, itmax
          xm    = 0.5_dp * (a+b)
          tol1  = tol*ABS(x)+zeps
          tol2  = 2.0_dp * tol1
          skip0 = .FALSE.
          skip1 = .FALSE.
          IF(ABS(x-xm).LE.(tol2-0.5_dp*(b-a))) EXIT
          IF(ABS(e).GT.tol1) THEN
             d1 = 2.0_dp*(b-a)
             d2 = d1
             IF(dw.NE.dx) d1=(w-x)*dx/(dx-dw)
             IF(dv.NE.dx) d2=(v-x)*dx/(dx-dv)
             u1   = x+d1
             u2   = x+d2
             ok1  = ((a-u1)*(u1-b).GT.0.0_dp).AND.(dx*d1.LE.0.0_dp)
             ok2  = ((a-u2)*(u2-b).GT.0.0_dp).AND.(dx*d2.LE.0.0_dp)
             olde = e
             e    = d
             IF(.NOT.(ok1.OR.ok2)) THEN
                skip0 = .TRUE.
             ELSE IF (ok1.AND.ok2) THEN
                IF(ABS(d1).LT.ABS(d2)) THEN
                   d = d1
                ELSE
                   d = d2
                ENDIF
             ELSE IF (ok1) THEN
                d = d1
             ELSE
                d = d2
             ENDIF
             IF (.NOT.skip0) THEN
                IF(ABS(d).GT.ABS(0.5_dp*olde)) skip0=.TRUE. 
                IF (.NOT.skip0) THEN
                   u=x+d
                   IF((u-a).LT.tol2 .OR. (b-u).LT.tol2) d=SIGN(tol1,xm-x)
                   skip1=.TRUE.
                END IF
             END IF
          ENDIF
          IF (.NOT.skip1) THEN
             IF(dx.GE.0.0_dp) THEN
                e=a-x
             ELSE
                e=b-x
             ENDIF
             d=0.5_dp*e
          END IF
          IF(ABS(d).GE.tol1) THEN
             u=x+d
             du=cg_deval1d(force_env,u,pcom,xicom,fu,error)
             loc_iter = loc_iter + 1
          ELSE
             u=x+SIGN(tol1,d)
             du=cg_deval1d(force_env,u,pcom,xicom,fu,error)
             loc_iter = loc_iter + 1
             IF(fu.GT.fx) EXIT
          ENDIF
          IF(fu.LE.fx) THEN
             IF(u.GE.x) THEN
                a=x
             ELSE
                b=x
             ENDIF
             v=w;   fv=fw;  dv=dw;  w=x
             fw=fx; dw=dx;  x=u;    fx=fu;  dx=du
          ELSE
             IF(u.LT.x) THEN
                a=u
             ELSE
                b=u
             ENDIF
             IF(fu.LE.fw .OR. w.EQ.x) THEN
                v=w; fv=fw; dv=dw
                w=u; fw=fu; dw=du
             ELSE IF(fu.LE.fv .OR. v.EQ.x .OR. v.EQ.w) THEN
                v=u
                fv=fu
                dv=du
             ENDIF
          ENDIF
       END DO
       IF (ionode) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
               "***","BRENT   - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
          IF (iter==itmax+1) &
               WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,T78,A)")&
               "***","BRENT - NUMBER OF ITERATIONS EXCEEDED ","***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       END IF
       CPPrecondition((iter/=itmax+1),cp_failure_level,routineP,error,failure)
       xmin=x
       dbrent=fx
    END IF
  END FUNCTION cg_dbrent
  
!!****f* cg_optimizer/cg_eval1d [1.0] *
!!
!!   NAME
!!     cg_eval1d
!!
!!   FUNCTION
!!     Evaluates energy in one dimensional space defined by the point
!!     pcom and with direction xicom, position x
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************  
  FUNCTION cg_eval1d(force_env,x,pcom,xicom,error) RESULT(my_val)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: my_val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cg_eval1d', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: stat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ALLOCATE(xvec(SIZE(pcom)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       xvec = pcom + x*xicom
       CALL get_energy_and_gradients(force_env, xvec, my_val)
       DEALLOCATE(xvec, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END FUNCTION cg_eval1d

!!****f* cg_optimizer/cg_deval1d [1.0] *
!!
!!   NAME
!!     cg_deval1d
!!
!!   FUNCTION
!!     Evaluates derivatives in one dimensional space defined by the point
!!     pcom and with direction xicom, position x
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino] 
!!
!!*** **********************************************************************
  FUNCTION cg_deval1d(force_env,x,pcom,xicom,fval,error) RESULT(my_val)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: my_val, fval
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cg_deval1d', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: stat
    REAL(KIND=dp)                            :: energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: grad, xvec
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ALLOCATE(xvec(SIZE(pcom)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(grad(SIZE(pcom)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       xvec = pcom + x * xicom
       CALL get_energy_and_gradients(force_env, xvec, energy, grad)
       my_val = DOT_PRODUCT(grad,xicom)
       fval   = energy
       DEALLOCATE(xvec, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(grad, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END FUNCTION cg_deval1d


!!****f* cg_optimizer/FindFit [1.0] *
!!
!!   NAME
!!     FindMin
!!
!!   FUNCTION
!!     Find the minimum of a parabolic function obtained with a least square fit
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [fawzi] 
!!
!!*** **********************************************************************
  FUNCTION FindMin(x,y,dy,error) RESULT(res)
    IMPLICIT NONE
    ! Arguments
    REAL(kind=dp), DIMENSION(:) :: x,y,dy
    REAL(kind=dp) :: res
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'findMin', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), DIMENSION(2*SIZE(x)) :: w,b
    REAL(kind=dp), DIMENSION(2*SIZE(x),3) ::f
    REAL(kind=dp) :: sum_x,sum_xx,spread,diag(3),u(2*SIZE(x),3),vt(3,3), res1(3),res2(3),tmpw(1),res3(3)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: work
    INTEGER :: i,np,min_pos,lwork, iwork(8*3), info, stat
    LOGICAL :: failure
    
    failure=.FALSE.
    np=SIZE(x)
    CPPrecondition(np>1,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       sum_x=0._dp
       sum_xx=0._dp
       min_pos=1
       DO i=1,np
          sum_xx=sum_xx+x(i)**2
          sum_x=sum_x+x(i)
          IF (y(min_pos)>y(i)) min_pos=i
       END DO
       spread=SQRT(sum_xx/REAL(np,dp)-(sum_x/REAL(np,dp))**2)
       DO i=1,np
          w(i)=EXP(-(REAL(np-i,dp))**2/(REAL(2*9,dp)))
          w(i+np)=2._dp*w(i)
       END DO
       DO i=1,np
          f(i,1)=w(i)
          f(i,2)=x(i)*w(i)
          f(i,3)=x(i)**2*w(i)
          f(i+np,1)=0
          f(i+np,2)=w(i+np)
          f(i+np,3)=2*x(i)*w(i+np)
       END DO
       DO i=1,np
          b(i)=y(i)*w(i)
          b(i+np)=dy(i)*w(i+np)
       END DO
       lwork=-1
       CALL dgesdd('S',SIZE(f,1),SIZE(f,2),f,SIZE(f,1),diag,u,SIZE(u,1),vt,SIZE(vt,1),tmpw,lwork,&
            iwork,info)
       lwork=ceiling(tmpw(1))
       ALLOCATE(work(lwork),stat=stat)
       CALL dgesdd('S',SIZE(f,1),SIZE(f,2),f,SIZE(f,1),diag,u,SIZE(u,1),vt,SIZE(vt,1),work,lwork,&
            iwork,info)
       DEALLOCATE(work, stat=stat)
       CALL dgemv('T',SIZE(u,1),SIZE(u,2),1._dp,u,SIZE(u,1),b,1,0._dp,res1,1)
       DO i=1,3
          res2(i)=res1(i)/diag(i)
       END DO
       CALL dgemv('T',3,3,1._dp,vt,SIZE(vt,1),res2,1,0._dp,res3,1)
       res=-0.5*res3(2)/res3(3)
    END IF
  END FUNCTION FindMin

END MODULE cg_optimizer 
