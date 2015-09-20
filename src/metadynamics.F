!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performs the metadynamics calculation
!> \par History
!>      01.2005 created [fawzi and ale]
!>      11.2007 Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE metadynamics
  USE bibliography,                    ONLY: VandenCic2006
#if defined (__PLUMED_CP2K) || defined (__PLUMED2)
  USE cell_types,                      ONLY: cell_type,&
                                             pbc_cp2k_plumed_getset_cell
#else
  USE cell_types,                      ONLY: cell_type
#endif
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             torsion_colvar_id
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: do_wall_m,&
                                             do_wall_p,&
                                             do_wall_reflective
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast, cp2k_is_parallel
  USE metadynamics_types,              ONLY: hills_env_type,&
                                             meta_env_type,&
                                             metavar_type,&
                                             multiple_walkers_type
  USE metadynamics_utils,              ONLY: add_hill_single,&
                                             get_meta_iter_level,&
                                             meta_walls,&
                                             restart_hills,&
                                             synchronize_multiple_walkers
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
#if defined (__PLUMED2)
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             femtoseconds,&
                                             joule,&
                                             kelvin,kjmol,&
                                             picoseconds
#else
  USE physcon,                         ONLY: boltzmann,&
                                             femtoseconds,&
                                             joule,&
                                             kelvin
#endif
  USE reference_manager,               ONLY: cite_reference
  USE simpar_types,                    ONLY: simpar_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics'

  PUBLIC :: metadyn_forces, metadyn_integrator
  PUBLIC :: metadyn_velocities_colvar, metadyn_write_colvar
  PUBLIC :: metadyn_initialise_plumed, metadyn_finalise_plumed

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param simpar ...
!> \param itimes ...
! *****************************************************************************
  SUBROUTINE metadyn_initialise_plumed(force_env, simpar, itimes)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar
    INTEGER, INTENT(IN)                      :: itimes

    CHARACTER(LEN=*), PARAMETER :: routineN = 'metadyn_initialise_plumed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

#if defined (__PLUMED_CP2K) || defined (__PLUMED2)
    INTEGER                                  :: natom_plumed
    REAL(KIND=dp)                            :: timestep_plumed
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
#endif

#if defined (__PLUMED_CP2K)
    INTEGER                                  :: nsp_plumed, i_kinds
    INTEGER, DIMENSION(90)                   :: na_plumed
    REAL(KIND=dp), DIMENSION(90)             :: pma0_plumed
#endif

#if defined (__PLUMED2)
    INTEGER :: plumedavailable
#endif

    CALL timeset(routineN,handle)
    CPASSERT(ASSOCIATED(simpar))
    CPASSERT(ASSOCIATED(force_env))

#if defined (__PLUMED_CP2K) || defined (__PLUMED2)
    NULLIFY(cell,subsys)
    CALL force_env_get(force_env, subsys=subsys, cell=cell)
    CALL pbc_cp2k_plumed_getset_cell (cell, set=.TRUE.) !Store the cell pointer for later use.
    timestep_plumed=simpar%dt
    natom_plumed=subsys%particles%n_els
#endif

#if defined (__PLUMED_CP2K)
    !in a parallel run plumed only needs to be initialized on the master node
    IF (force_env%para_env%ionode .EQV. .TRUE.) THEN
        nsp_plumed=subsys%atomic_kinds%n_els
        DO i_kinds=1,nsp_plumed
            na_plumed(i_kinds) = subsys%atomic_kinds%els(i_kinds)%natom
            pma0_plumed(i_kinds) = subsys%atomic_kinds%els(i_kinds)%mass
        END DO
        CALL INIT_METADYN(natom_plumed, nsp_plumed, na_plumed, &
            timestep_plumed, itimes, pma0_plumed, force_env%meta_env%plumed_input_file)
    END IF
#endif

#if defined (__PLUMED2)
    CALL plumed_f_installed(plumedavailable)

    IF(plumedavailable<=0)THEN
        WRITE(*,*)"NO PLUMED IN YOUR LD_LIBRARY_PATH?"
        STOP
    ELSE
        CALL plumed_f_gcreate()
        IF (cp2k_is_parallel) CALL plumed_f_gcmd("setMPIFComm"//CHAR(0),force_env%para_env%group)
        CALL plumed_f_gcmd("setRealPrecision"//CHAR(0),8)
        CALL plumed_f_gcmd("setMDEnergyUnits"//CHAR(0),kjmol) ! Hartree to kjoule/mol
        CALL plumed_f_gcmd("setMDLengthUnits"//CHAR(0),angstrom*0.1_dp) ! Bohr to nm
        CALL plumed_f_gcmd("setMDTimeUnits"//CHAR(0),picoseconds) !  atomic units to ps
        CALL plumed_f_gcmd("setPlumedDat"//CHAR(0),force_env%meta_env%plumed_input_file//CHAR(0))
        CALL plumed_f_gcmd("setLogFile"//CHAR(0),"PLUMED.OUT"//CHAR(0))
        CALL plumed_f_gcmd("setNatoms"//CHAR(0),natom_plumed)
        CALL plumed_f_gcmd("setMDEngine"//CHAR(0),"cp2k"//CHAR(0))
        CALL plumed_f_gcmd("setTimestep"//CHAR(0),timestep_plumed)
        CALL plumed_f_gcmd("init"//CHAR(0),0)
    END IF
#endif

#if !defined (__PLUMED_CP2K) && !defined (__PLUMED2)
    CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"Requested to use plumed for metadynamics, but cp2k was"//&
                   " not compiled with plumed support."//&
                   CPSourceFileRef,&
                   failure)
#endif

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief makes sure PLUMED is shut down cleanly 
! *****************************************************************************
SUBROUTINE metadyn_finalise_plumed()

    CHARACTER(LEN=*), PARAMETER :: routineN = 'metadyn_finalise_plumed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

#if defined (__PLUMED2)
    CALL plumed_f_gfinalize()
#endif

    CALL timestop(handle)

END SUBROUTINE

! *****************************************************************************
!> \brief  General driver for applying metadynamics
!> \param force_env ...
!> \param itimes ...
!> \param vel ...
!> \param rand ...
!> \date   01.2009
!> \par History
!>      01.2009 created
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE metadyn_integrator(force_env, itimes, vel, rand)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: vel
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: rand

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_integrator', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, plumed_itimes
    LOGICAL                                  :: failure
#if defined (__PLUMED_CP2K) || defined (__PLUMED2)
    INTEGER                                  :: i_part, natom_plumed
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
#endif
#if defined (__PLUMED_CP2K)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: force_plumed, pos_plumed
#endif
#if defined (__PLUMED2)
    TYPE(meta_env_type), POINTER             :: meta_env
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mass_plumed
    REAL(KIND=dp), DIMENSION(3,3)            :: plu_vir, celltbt
    REAL(KIND=dp)                            :: stpcfg
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: pos_plumed_x, pos_plumed_y, pos_plumed_z
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: force_plumed_x, force_plumed_y, force_plumed_z
#endif

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ! Apply Metadynamics
    IF (ASSOCIATED(force_env%meta_env)) THEN
        IF (force_env%meta_env%use_plumed .EQV. .TRUE.) THEN
            plumed_itimes = itimes
#if defined (__PLUMED_CP2K)
            CALL force_env_get(force_env, subsys=subsys,cell=cell)
            natom_plumed=subsys%particles%n_els
            ALLOCATE (pos_plumed(1:3,natom_plumed))
            ALLOCATE (force_plumed(1:3,natom_plumed))
            force_plumed(:,:)=0.0d0

            DO i_part=1,natom_plumed
              pos_plumed(:,i_part)=subsys%particles%els(i_part)%r(:)
            END DO
           
            !in a parallel run plumed colvars are only processed on the
            !master node and then broadcast to the others.
            IF (force_env%para_env%ionode .EQV. .TRUE.) THEN
                CALL META_FORCE_CALCULATION(pos_plumed,force_plumed,natom_plumed,cell)
            END IF
            CALL mp_bcast (force_plumed,                     &
                           force_env%para_env%source,        &
                           force_env%para_env%group)
            
            DO i_part=1,natom_plumed
              subsys%particles%els(i_part)%f(:) = subsys%particles%els(i_part)%f(:) + force_plumed(:,i_part)
            END DO

            DEALLOCATE (force_plumed, pos_plumed)
#endif
#if defined (__PLUMED2)
            CALL force_env_get(force_env, meta_env=meta_env, subsys=subsys,cell=cell)
            natom_plumed=subsys%particles%n_els
            ALLOCATE (pos_plumed_x(natom_plumed))
            ALLOCATE (pos_plumed_y(natom_plumed))
            ALLOCATE (pos_plumed_z(natom_plumed))
        
            ALLOCATE (force_plumed_x(natom_plumed))
            ALLOCATE (force_plumed_y(natom_plumed))
            ALLOCATE (force_plumed_z(natom_plumed))

            ALLOCATE (mass_plumed(natom_plumed))

            force_plumed_x(:)=0.0d0 
            force_plumed_y(:)=0.0d0
            force_plumed_z(:)=0.0d0

            DO i_part=1,natom_plumed
              pos_plumed_x(i_part)=subsys%particles%els(i_part)%r(1)
              pos_plumed_y(i_part)=subsys%particles%els(i_part)%r(2)
              pos_plumed_z(i_part)=subsys%particles%els(i_part)%r(3)
              mass_plumed(i_part)=subsys%particles%els(i_part)%atomic_kind%mass
            END DO

            ! stupid cell conversion, to be fixed
            celltbt(1,1)=cell%hmat(1,1)
            celltbt(1,2)=cell%hmat(1,2)
            celltbt(1,3)=cell%hmat(1,3)
            celltbt(2,1)=cell%hmat(2,1)
            celltbt(2,2)=cell%hmat(2,2)
            celltbt(2,3)=cell%hmat(2,3)
            celltbt(3,1)=cell%hmat(3,1)
            celltbt(3,2)=cell%hmat(3,2)
            celltbt(3,3)=cell%hmat(3,3)
  
            ! virial
            plu_vir(:,:)=0.0d0
  
            CALL force_env_get(force_env, potential_energy=stpcfg)
  
            CALL plumed_f_gcmd("setStep"//CHAR(0),plumed_itimes)      
            CALL plumed_f_gcmd("setPositionsX"//CHAR(0),pos_plumed_x(:)) 
            CALL plumed_f_gcmd("setPositionsY"//CHAR(0),pos_plumed_y(:))
            CALL plumed_f_gcmd("setPositionsZ"//CHAR(0),pos_plumed_z(:))
            CALL plumed_f_gcmd("setMasses"//CHAR(0), mass_plumed(:)) 
            CALL plumed_f_gcmd("setBox"//CHAR(0),celltbt)
            CALL plumed_f_gcmd("setVirial"//CHAR(0),plu_vir)    ! PluMed Output, NOT ACTUALLY USED !
            CALL plumed_f_gcmd("setEnergy"//CHAR(0),stpcfg) 
            CALL plumed_f_gcmd("setForcesX"//CHAR(0),force_plumed_x(:)) ! PluMed Output !
            CALL plumed_f_gcmd("setForcesY"//CHAR(0),force_plumed_y(:)) ! PluMed Output !
            CALL plumed_f_gcmd("setForcesZ"//CHAR(0),force_plumed_z(:)) ! PluMed Output !
            
            CALL plumed_f_gcmd("prepareCalc"//CHAR(0),0) 
            CALL plumed_f_gcmd("prepareDependencies"//CHAR(0),0)
            CALL plumed_f_gcmd("shareData"//CHAR(0),0)
            
            CALL plumed_f_gcmd("performCalc"//CHAR(0),0)

            DO i_part=1,natom_plumed
               subsys%particles%els(i_part)%f(1) = subsys%particles%els(i_part)%f(1) + force_plumed_x(i_part)
               subsys%particles%els(i_part)%f(2) = subsys%particles%els(i_part)%f(2) + force_plumed_y(i_part)
               subsys%particles%els(i_part)%f(3) = subsys%particles%els(i_part)%f(3) + force_plumed_z(i_part)
            END DO
 
            DEALLOCATE (pos_plumed_x, pos_plumed_y, pos_plumed_z)
            DEALLOCATE (force_plumed_x, force_plumed_y, force_plumed_z)
            DEALLOCATE (mass_plumed)

            ! Constraints ONLY of Fixed Atom type
            CALL fix_atom_control(force_env)
#endif

#if !defined (__PLUMED_CP2K) && !defined (__PLUMED2)
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"Requested to use plumed for metadynamics, but cp2k was"//&
                   " not compiled with plumed support."//&
                   CPSourceFileRef,&
                   failure)
#endif


        ELSE
            IF(force_env%meta_env%langevin) THEN
              IF (.NOT.PRESENT(rand)) THEN
                 CALL cp_unimplemented_error(fromWhere=routineP, &
                      message="Langevin on COLVAR not implemented for this MD ensemble!", &
                      error_level=cp_fatal_level)
              END IF
              !    *** Velocity Verlet for Langevin S(t)->S(t+1)
              CALL metadyn_position_colvar(force_env)
              !    *** Forces from Vs and  S(X(t+1))
              CALL metadyn_forces(force_env)
              !    *** Velocity Verlet for Langeving *** v(t+1/2)--> v(t)
              CALL metadyn_velocities_colvar(force_env, rand)
            ELSE
              CALL metadyn_forces(force_env, vel)
            ENDIF
            !    *** Write down COVAR informations
            CALL metadyn_write_colvar(force_env)
        END IF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE metadyn_integrator

! *****************************************************************************
!> \brief add forces to the subsys due to the metadynamics run
!>      possibly modifies the velocites (if reflective walls are applied)
!> \param force_env ...
!> \param vel ...
!> \par History
!>      04.2004 created
! *****************************************************************************
  SUBROUTINE metadyn_forces(force_env,vel)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: vel

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, icolvar, ii, &
                                                iwall
    LOGICAL                                  :: explicit, failure
    REAL(kind=dp)                            :: check_val, diff_ss, dt, &
                                                ekin_w, fac_t, fft, norm, &
                                                rval, scal, scalf, ss0_test, &
                                                tol_ekin
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ss0_section, vvp_section

    failure=.FALSE.
    NULLIFY(logger,meta_env)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    NULLIFY(colvar_p,subsys,cv,ss0_section, vvp_section)
    CALL force_env_get(force_env, subsys=subsys)

    dt = meta_env%dt
    IF (.NOT.meta_env%restart) meta_env%n_steps=meta_env%n_steps+1

    ! Initialize velocity
    IF (meta_env%restart.AND.meta_env%extended_lagrange) THEN
       meta_env%ekin_s = 0.0_dp
       DO i_c=1,meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%vvp = next_random_number(force_env%globenv%gaussian_rng_stream)
          meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       END DO
       ekin_w = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
       fac_t  = SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
       DO i_c = 1, meta_env%n_colvar
          cv     => meta_env%metavar(i_c)
          cv%vvp = cv%vvp*fac_t
       ENDDO
       meta_env%ekin_s = 0.0_dp
    END IF

          !    *** Velocity Verlet for Langevin S(t)->S(t+1)
    ! compute ss and the derivative of ss with respect to the atomic positions
    DO i_c=1,meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       icolvar = cv%icolvar
       CALL colvar_eval_glob_f(icolvar,force_env)
       cv%ss       = subsys%colvar_p(icolvar)%colvar%ss

       ! Setup the periodic flag if the COLVAR is (-pi,pi] periodic
       cv%periodic = (subsys%colvar_p(icolvar)%colvar%type_id==torsion_colvar_id)

       ! Restart for Extended Lagrangian Metadynamics
       IF (meta_env%restart) THEN
          ! Initialize the position of the collective variable in the extended lagrange
          ss0_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS0")
          CALL section_vals_get(ss0_section, explicit=explicit)
          IF (explicit) THEN
             CALL section_vals_val_get(ss0_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval)
             cv%ss0 = rval
          ELSE
             cv%ss0 = cv%ss
          END IF
          vvp_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_VVP")
          CALL section_vals_get(vvp_section, explicit=explicit)
          IF (explicit) THEN
             CALL section_vals_val_get(vvp_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval)
             cv%vvp = rval
          END IF
       END IF
       !
       IF(.NOT.meta_env%extended_lagrange) THEN
          cv%ss0 = cv%ss
          cv%vvp = 0.0_dp
       END IF
    ENDDO
    ! History dependent forces (evaluated at s0)
    IF(meta_env%do_hills) CALL  hills(meta_env)

    ! Apply walls to the colvars
    CALL meta_walls(meta_env)

    meta_env%restart = .FALSE.
    IF(.NOT.meta_env%extended_lagrange)THEN
       meta_env%ekin_s = 0.0_dp
       meta_env%epot_s = 0.0_dp
       meta_env%epot_walls = 0.0_dp
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%epot_s = 0.0_dp
          cv%ff_s   = 0.0_dp
          meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
          icolvar=cv%icolvar
          NULLIFY(particles)
          CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
               particles=particles)
          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             fft = cv%ff_hills + cv%ff_walls
             particles%els(i)%f=particles%els(i)%f+fft*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO
       ENDDO
    ELSE
       meta_env%ekin_s = 0.0_dp
       meta_env%epot_s = 0.0_dp
       meta_env%epot_walls = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          diff_ss = cv%ss-cv%ss0
          IF (cv%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          cv%epot_s = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
          cv%ff_s   = cv%lambda*(diff_ss)
          icolvar   = cv%icolvar
          ! forces on the atoms
          NULLIFY(particles)
          CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
               particles=particles)
          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f- cv%ff_s*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO
          !  velocity verlet on the s0 if NOT langevin
          IF (.NOT.meta_env%langevin) THEN
             fft    = cv%ff_s + cv%ff_hills + cv%ff_walls
             cv%vvp = cv%vvp+dt*fft/cv%mass
             meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
             meta_env%epot_s = meta_env%epot_s + cv%epot_s
             meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
          END IF
       ENDDO
       !  velocity rescaling on the s0
       IF(meta_env%tempcontrol.AND.(.NOT.meta_env%langevin))THEN
          ekin_w   = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
          tol_ekin = 0.5_dp*meta_env%toll_temp*REAL(meta_env%n_colvar,KIND=dp)
          IF(ABS(ekin_w-meta_env%ekin_s)> tol_ekin)THEN
             fac_t=SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
             DO i_c = 1, meta_env%n_colvar
                cv     => meta_env%metavar(i_c)
                cv%vvp = cv%vvp*fac_t
             ENDDO
             meta_env%ekin_s = ekin_w
          ENDIF
       ENDIF
       ! Reflective Wall only for s0
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          IF (cv%do_wall) THEN
             DO iwall = 1, SIZE(cv%walls)
                SELECT CASE(cv%walls(iwall)%id_type)
                CASE(do_wall_reflective)
                   ss0_test=cv%ss0+dt*cv%vvp
                   IF (cv%periodic) THEN
                      ! A periodic COLVAR is always within [-pi,pi]
                      ss0_test = SIGN(1.0_dp,ASIN(SIN(ss0_test)))*ACOS(COS(ss0_test))
                   END IF
                   SELECT CASE(cv%walls(iwall)%id_direction)
                   CASE(do_wall_p)
                      IF ((ss0_test > cv%walls(iwall)%pos).AND.(cv%vvp>0)) cv%vvp = -cv%vvp
                   CASE(do_wall_m)
                      IF ((ss0_test < cv%walls(iwall)%pos).AND.(cv%vvp<0)) cv%vvp = -cv%vvp
                   END SELECT
                END SELECT
             END DO
          ENDIF
       ENDDO
       ! Update of ss0 if NOT langevin
       IF (.NOT.meta_env%langevin) THEN
          DO i_c = 1, meta_env%n_colvar
             cv => meta_env%metavar(i_c)
             cv%ss0 = cv%ss0 + dt*cv%vvp
             IF (cv%periodic) THEN
                ! A periodic COLVAR is always within [-pi,pi]
                cv%ss0 = SIGN(1.0_dp,ASIN(SIN(cv%ss0)))*ACOS(COS(cv%ss0))
             END IF
          ENDDO
       END IF
    ENDIF
    ! Constraints ONLY of Fixed Atom type
    CALL fix_atom_control(force_env)

    ! Reflective Wall only for ss
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       IF (cv%do_wall) THEN
          DO iwall = 1, SIZE(cv%walls)
             SELECT CASE(cv%walls(iwall)%id_type)
             CASE(do_wall_reflective)
                SELECT CASE(cv%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (cv%ss < cv%walls(iwall)%pos) CYCLE
                   check_val = -1.0_dp
                CASE(do_wall_m)
                   IF (cv%ss > cv%walls(iwall)%pos) CYCLE
                   check_val =  1.0_dp
                END SELECT
                NULLIFY(particles)
                icolvar=cv%icolvar
                CALL cp_subsys_get(subsys, colvar_p=colvar_p, particles=particles)
                scal=0.0_dp
                scalf=0.0_dp
                norm=0.0_dp
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      scal=scal+vel(1,i)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                      scal=scal+vel(2,i)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                      scal=scal+vel(3,i)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   ELSE
                      scal=scal+particles%els(i)%v(1)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                      scal=scal+particles%els(i)%v(2)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                      scal=scal+particles%els(i)%v(3)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   END IF
                   scalf=scalf+particles%els(i)%f(1)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                   scalf=scalf+particles%els(i)%f(2)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                   scalf=scalf+particles%els(i)%f(3)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(1,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(2,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(3,ii)**2
                ENDDO
                IF (norm /= 0.0_dp) scal = scal/norm
                IF (norm /= 0.0_dp) scalf= scalf/norm

                IF (scal*check_val>0.0_dp) CYCLE
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      vel(:,i)=vel(:,i)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   ELSE
                      particles%els(i)%v(:)=particles%els(i)%v(:)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   END IF
                   ! Nullify forces along the colvar (this avoids the weird behaviors of the reflective wall)
                    particles%els(i)%f(:)=particles%els(i)%f(:)-colvar_p(icolvar)%colvar%dsdr(:,ii)*scalf
                ENDDO
             END SELECT
          END DO
       END IF
    ENDDO

    CALL timestop(handle)
  END SUBROUTINE metadyn_forces

! *****************************************************************************
!> \brief Evolves velocities COLVAR according to
!>        Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \param force_env ...
!> \param rand ...
!> \date  01.2009
!> \author Fabio Sterpone and Teodoro Laino
! *****************************************************************************
  SUBROUTINE metadyn_velocities_colvar(force_env,rand)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: rand

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_velocities_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: diff_ss, dt, fft, sigma
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    ! Add citation
    IF (meta_env%langevin) CALL cite_reference(VandenCic2006)

    dt = meta_env%dt
    ! History dependent forces (evaluated at s0)
    IF(meta_env%do_hills) CALL  hills(meta_env)

    ! Evolve Velocities
    meta_env%ekin_s = 0.0_dp
    meta_env%epot_walls = 0.0_dp
    DO i_c= 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       diff_ss = cv%ss-cv%ss0
       IF (cv%periodic) THEN
          ! The difference of a periodic COLVAR is always within [-pi,pi]
          diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
       END IF
       cv%epot_s = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
       cv%ff_s   = cv%lambda*(diff_ss)

       fft    = cv%ff_s+cv%ff_hills
       sigma  = SQRT((meta_env%temp_wanted*kelvin)*2.0_dp*(boltzmann/joule)*cv%gamma/cv%mass)
       cv%vvp = cv%vvp+0.5_dp*dt*fft/cv%mass-0.5_dp*dt*cv%gamma*cv%vvp+&
                0.5_dp*SQRT(dt)*sigma*rand(i_c)
       meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
    ENDDO
    CALL timestop(handle)

  END SUBROUTINE metadyn_velocities_colvar

! *****************************************************************************
!> \brief Evolves COLVAR position
!>        Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \param force_env ...
!> \date  01.2009
!> \author Fabio Sterpone and Teodoro Laino
! *****************************************************************************
  SUBROUTINE metadyn_position_colvar(force_env)
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_position_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()

    ! Add citation
    IF (meta_env%langevin) CALL cite_reference(VandenCic2006)
    dt = meta_env%dt

    ! Update of ss0
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       cv%ss0 =cv%ss0 +dt*cv%vvp
       IF (cv%periodic) THEN
          ! A periodic COLVAR is always within [-pi,pi]
          cv%ss0 = SIGN(1.0_dp,ASIN(SIN(cv%ss0)))*ACOS(COS(cv%ss0))
       END IF
    ENDDO
    CALL timestop(handle)

  END SUBROUTINE metadyn_position_colvar

! *****************************************************************************
!> \brief Write down COLVAR information evolved  according to
!>        Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \param force_env ...
!> \date  01.2009
!> \author Fabio Sterpone and Teodoro Laino
! *****************************************************************************
  SUBROUTINE metadyn_write_colvar(force_env)
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_write_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: diff_ss, temp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()

    ! If Langevin we need to recompute few quantities
    ! This does not apply to the standard lagrangian scheme since it is
    ! implemented with a plain Newton integration scheme.. while Langevin
    ! follows the correct Verlet integration.. This will have to be made
    ! uniform in the future (Teodoro Laino - 01.2009)
    IF (meta_env%langevin) THEN
       meta_env%ekin_s = 0.0_dp
       meta_env%epot_s = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          diff_ss = cv%ss-cv%ss0
          IF (cv%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          cv%epot_s = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
          cv%ff_s   = cv%lambda*(diff_ss)

          meta_env%epot_s = meta_env%epot_s + cv%epot_s
          meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       ENDDO
    END IF

    ! write COLVAR file
    iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
         "PRINT%COLVAR",extension=".metadynLog")
    IF (iw>0) THEN
       IF (meta_env%extended_lagrange) THEN
          WRITE(iw,'(f16.8,70f15.8)')meta_env%time*femtoseconds, &
               (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ss,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_s,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_walls,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%vvp,i=1,meta_env%n_colvar), &
                meta_env%epot_s, &
                meta_env%hills_env%energy, &
                meta_env%epot_walls, &
               (meta_env%ekin_s)*2.0_dp/(REAL(meta_env%n_colvar,KIND=dp))*kelvin
       ELSE
          WRITE(iw,'(f16.8,40f13.5)') meta_env%time*femtoseconds,&
               (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar),&
               (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar),&
               (meta_env%metavar(i)%ff_walls,i=1,meta_env%n_colvar),&
                meta_env%hills_env%energy,&
                meta_env%epot_walls
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
         "PRINT%COLVAR")

    ! Temperature for COLVAR
    IF (meta_env%extended_lagrange) THEN
       temp = meta_env%ekin_s*2.0_dp/(REAL(meta_env%n_colvar,KIND=dp))*kelvin
       meta_env%avg_temp = (meta_env%avg_temp*REAL(meta_env%n_steps,KIND=dp)+&
            temp)/REAL(meta_env%n_steps+1,KIND=dp)
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%TEMPERATURE_COLVAR",extension=".metadynLog")
       IF (iw > 0) THEN
          WRITE (iw, '(T2,79("-"))')
          WRITE (iw,'( A,T51,f10.2,T71,f10.2)' )' COLVARS INSTANTANEOUS/AVERAGE TEMPERATURE ',&
               temp, meta_env%avg_temp
          WRITE (iw, '(T2,79("-"))')
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%TEMPERATURE_COLVAR")
    END IF
    CALL timestop(handle)

  END SUBROUTINE metadyn_write_colvar

! *****************************************************************************
!> \brief Major driver for adding hills and computing forces due to the history
!>        dependent term
!> \param meta_env ...
!> \par History
!>      04.2004 created
!>      10.2008 Teodoro Laino [tlaino] - University of Zurich
!>              Major rewriting and addition of multiple walkers
! *****************************************************************************
  SUBROUTINE hills(meta_env)
    TYPE(meta_env_type), POINTER             :: meta_env

    CHARACTER(len=*), PARAMETER :: routineN = 'hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_hills, ih, &
                                                intermeta_steps, iter_nr, iw, &
                                                n_colvar, n_hills_start, &
                                                n_step, stat
    LOGICAL                                  :: failure, force_gauss
    REAL(KIND=dp)                            :: cut_func, dfunc, diff, dp2, &
                                                frac, slow_growth, &
                                                V_now_here, V_to_fes, wtww, ww
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ddp, denf, diff_ss, &
                                                local_last_hills, numf
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(multiple_walkers_type), POINTER     :: multiple_walkers

    CALL timeset(routineN,handle)

    failure   = .FALSE.
    NULLIFY(hills_env, multiple_walkers, logger, colvars, ddp, local_last_hills)
    hills_env => meta_env%hills_env
    logger    => cp_get_default_logger()
    colvars   => meta_env%metavar
    n_colvar  =  meta_env%n_colvar
    n_step    =  meta_env%n_steps

    ! Create a temporary logger level specific for metadynamics
    CALL cp_add_iter_level(logger%iter_info,"METADYNAMICS")
    CALL get_meta_iter_level(meta_env, iter_nr)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_nr)

    ! Set-up restart if any
    IF (meta_env%hills_env%restart) THEN
       meta_env%hills_env%restart = .FALSE.
       IF (meta_env%well_tempered) THEN
          CALL restart_hills(hills_env%ss_history, hills_env%delta_s_history, hills_env%ww_history,&
               hills_env%ww, hills_env%n_hills, n_colvar, colvars, meta_env%metadyn_section,&
               invdt_history=hills_env%invdt_history)
       ELSE
          CALL restart_hills(hills_env%ss_history, hills_env%delta_s_history, hills_env%ww_history,&
               hills_env%ww, hills_env%n_hills, n_colvar, colvars, meta_env%metadyn_section)
       END IF
    END IF


    ! Proceed with normal calculation
    intermeta_steps=n_step-hills_env%old_hill_step
    force_gauss=.FALSE.
    IF ((hills_env%min_disp > 0.0_dp).AND.(hills_env%old_hill_number > 0).AND.&
        (intermeta_steps >= hills_env%min_nt_hills)) THEN
       ALLOCATE (ddp(meta_env%n_colvar), stat=stat)
       CPASSERT(stat==0)
       ALLOCATE (local_last_hills(meta_env%n_colvar), stat=stat)
       CPASSERT(stat==0)

       local_last_hills(1:n_colvar) = hills_env%ss_history(1:n_colvar,hills_env%old_hill_number)

       !RG Calculate the displacement
       dp2 = 0.0_dp
       DO i = 1, n_colvar
          ddp(i) = colvars(i)%ss0 - local_last_hills(i)
          IF (colvars(i)%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             ddp(i) = SIGN(1.0_dp,ASIN(SIN(ddp(i))))*ACOS(COS(ddp(i)))
          END IF
          dp2 = dp2 + ddp(i)**2
       ENDDO
       dp2 = SQRT(dp2)
       IF (dp2 > hills_env%min_disp) THEN
          force_gauss = .TRUE.
          iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
                                    "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog")
          IF (iw > 0) THEN
             WRITE (UNIT=iw,FMT="(A,F0.6,A,F0.6)")&
              " METADYN| Threshold value for COLVAR displacement exceeded: ",&
              dp2," > ",hills_env%min_disp
          END IF
          CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
                                            "PRINT%PROGRAM_RUN_INFO")
       END IF
       DEALLOCATE (ddp, stat=stat)
       CPASSERT(stat==0)
       DEALLOCATE (local_last_hills, stat=stat)
       CPASSERT(stat==0)
    END IF

    !RG keep into account adaptive hills
    IF (((MODULO(intermeta_steps,hills_env%nt_hills)==0).OR.force_gauss)&
        .AND.(.NOT.meta_env%restart).AND.(hills_env%nt_hills>0)) THEN
       IF (meta_env%do_multiple_walkers) multiple_walkers => meta_env%multiple_walkers

       n_hills_start = hills_env%n_hills
       ! Add the hill corresponding to this location
       IF (meta_env%well_tempered) THEN
          ! Well-Tempered scaling of hills height
          V_now_here = 0._dp
          DO ih=1,hills_env%n_hills
             dp2 = 0._dp
             DO i = 1, n_colvar
               diff = colvars(i)%ss0-hills_env%ss_history(i,ih)
               IF (colvars(i)%periodic) THEN
                  ! The difference of a periodic COLVAR is always within [-pi,pi]
                  diff = SIGN(1.0_dp,ASIN(SIN(diff)))*ACOS(COS(diff))
               END IF
               diff = (diff)/hills_env%delta_s_history(i,ih)
               dp2 = dp2 + diff**2
             ENDDO
             V_to_fes   = 1.0_dp + meta_env%wttemperature*hills_env%invdt_history(ih)
             V_now_here = V_now_here+hills_env%ww_history(ih)/V_to_fes*EXP(-0.5_dp*dp2)
          ENDDO
          wtww=hills_env%ww*EXP(-V_now_here*meta_env%invdt)
          ww=wtww*(1.0_dp + meta_env%wttemperature*meta_env%invdt)
          CALL add_hill_single(hills_env, colvars, ww, hills_env%n_hills, n_colvar,meta_env%invdt)
       ELSE
          CALL add_hill_single(hills_env, colvars, hills_env%ww, hills_env%n_hills, n_colvar)
       END IF
       ! Update local n_hills counter
       IF (meta_env%do_multiple_walkers) multiple_walkers%n_hills_local = multiple_walkers%n_hills_local+1

       hills_env%old_hill_number=hills_env%n_hills
       hills_env%old_hill_step=n_step

       ! Update iteration level for printing
       CALL get_meta_iter_level(meta_env, iter_nr)
       CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_nr)

       ! Print just program_run_info
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog")
       IF (iw>0) THEN
          IF (meta_env%do_multiple_walkers) THEN
             WRITE(iw,'(/,1X,"METADYN|",A,I0,A,I0,A,/)')&
                  ' Global/Local Hills number (',hills_env%n_hills,'/',multiple_walkers%n_hills_local,&
                  ') added.'
          ELSE
             WRITE(iw,'(/,1X,"METADYN|",A,I0,A,/)')' Hills number (',hills_env%n_hills,') added.'
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO")

       ! Handle Multiple Walkers
       IF (meta_env%do_multiple_walkers) THEN
          ! Print Local Hills file if requested
          iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
               "PRINT%HILLS",middle_name="LOCAL",extension=".metadynLog")
          IF (iw>0) THEN
             WRITE(iw,'(f12.1,30f13.5)')meta_env%time*femtoseconds,&
                  (hills_env%ss_history(ih,hills_env%n_hills),ih=1,n_colvar),&
                  (hills_env%delta_s_history(ih,hills_env%n_hills),ih=1,n_colvar),&
                   hills_env%ww_history(hills_env%n_hills)
          END IF
          CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
               "PRINT%HILLS")

          ! Check the communication buffer of the other walkers
          CALL synchronize_multiple_walkers(multiple_walkers, hills_env, colvars,&
               n_colvar, meta_env%metadyn_section)
       END IF

       ! Print Hills file if requested (for multiple walkers this file includes
       ! the Hills coming from all the walkers).
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%HILLS",extension=".metadynLog")
       IF (iw>0) THEN
          DO i_hills = n_hills_start+1, hills_env%n_hills
             WRITE(iw,'(f12.1,30f13.5)')meta_env%time*femtoseconds,&
                  (hills_env%ss_history(ih,i_hills),ih=1,n_colvar),&
                  (hills_env%delta_s_history(ih,i_hills),ih=1,n_colvar),&
                   hills_env%ww_history(i_hills)
          END DO
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%HILLS")
    END IF

    ! Computes forces due to the hills: history dependent term
    ALLOCATE(ddp(meta_env%n_colvar), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(diff_ss(meta_env%n_colvar), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(numf(meta_env%n_colvar), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(denf(meta_env%n_colvar), stat=stat)
    CPASSERT(stat==0)

    hills_env%energy = 0.0_dp
    DO ih=1,n_colvar
       colvars(ih)%ff_hills = 0.0_dp
    ENDDO
    DO ih=1,hills_env%n_hills
       slow_growth = 1.0_dp
       IF(hills_env%slow_growth .AND. (ih == hills_env%n_hills)) THEN
        slow_growth = REAL(n_step-hills_env%old_hill_step,dp)/REAL(hills_env%nt_hills,dp)
       END IF
       dp2 = 0._dp
       cut_func = 1.0_dp
       DO i = 1, n_colvar
          diff_ss(i) = colvars(i)%ss0-hills_env%ss_history(i,ih)
          IF (colvars(i)%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss(i) = SIGN(1.0_dp,ASIN(SIN(diff_ss(i))))*ACOS(COS(diff_ss(i)))
          END IF
          IF (hills_env%delta_s_history(i,ih) == 0.0_dp) THEN
            ! trick: scale = 0 is interpreted as infinitely wide Gaussian hill
            ! instead of infinitely narrow. This way one can combine several
            ! one-dimensional bias potentials in a multi-dimensional metadyn
            ! simulation.
            ddp(i) = 0.0_dp
          ELSE
            ddp(i) = (diff_ss(i))/hills_env%delta_s_history(i,ih)
          END IF
          dp2 = dp2 + ddp(i)**2
          IF(hills_env%tail_cutoff>0.0_dp) THEN
             frac = ABS(ddp(i))/hills_env%tail_cutoff
             numf(i) =  frac**hills_env%p_exp
             denf(i) =  frac**hills_env%q_exp
             cut_func =  cut_func * (1.0_dp-numf(i))/(1.0_dp-denf(i))
          ENDIF
       ENDDO
       ! ff_hills contains the "force" due to the hills
       dfunc = hills_env%ww_history(ih)*EXP(-0.5_dp*dp2)*slow_growth
       IF (meta_env%well_tempered) dfunc=dfunc/(1.0_dp + meta_env%wttemperature*hills_env%invdt_history(ih))
       hills_env%energy = hills_env%energy + dfunc * cut_func
       DO i=1,n_colvar
          IF (hills_env%delta_s_history(i,ih) /= 0.0_dp) THEN
            ! only apply a force when the Gaussian hill has a finite width in
            ! this direction
            colvars(i)%ff_hills = colvars(i)%ff_hills +&
                 ddp(i)/hills_env%delta_s_history(i,ih)*dfunc*cut_func
             IF(hills_env%tail_cutoff>0.0_dp .AND. ABS(diff_ss(i)) > 10.E-5_dp) THEN
               colvars(i)%ff_hills = colvars(i)%ff_hills +&
                 (hills_env%p_exp*numf(i)/(1.0_dp-numf(i))-hills_env%q_exp*denf(i)/(1.0_dp-denf(i)))*&
                  dfunc*cut_func/ABS(diff_ss(i))
             END IF
          END IF
       ENDDO
    ENDDO
    DEALLOCATE(ddp, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(diff_ss, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(numf, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(denf, stat=stat)
    CPASSERT(stat==0)

    CALL cp_rm_iter_level(logger%iter_info,"METADYNAMICS")

    CALL timestop(handle)

  END SUBROUTINE hills

END MODULE metadynamics
