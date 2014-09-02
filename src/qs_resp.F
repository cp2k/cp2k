!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief provides a resp fit for gas phase systems
!> \par History
!>      created
!>      Dorothea Golze [06.2012] (1) extension to periodic systems
!>                               (2) re-structured the code
!> \author Joost VandeVondele (02.2007)
! *****************************************************************************
MODULE qs_resp
  USE atomic_charges,                  ONLY: print_atomic_charges
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                             cp_unit_to_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: do_resp_minus_x_dir,&
                                             do_resp_minus_y_dir,&
                                             do_resp_minus_z_dir,&
                                             do_resp_x_dir,&
                                             do_resp_y_dir,&
                                             do_resp_z_dir,&
                                             use_perd_none,&
                                             use_perd_xyz
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_get_lval,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_wait
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE qs_collocate_density,            ONLY: calculate_rho_resp_all,&
                                             calculate_rho_resp_single
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_resp'

  PUBLIC :: resp_fit

  TYPE resp_type
   LOGICAL                                :: equal_charges, itc,& 
                                             nonperiodic_sys, rheavies
   INTEGER                                :: nres, ncons,&
                                             nrest_sec, ncons_sec,&
                                             npoints, stride(3), my_fit,&
                                             npoints_proc
   INTEGER, DIMENSION(:), POINTER         :: atom_surf_list
   INTEGER, DIMENSION(:,:), POINTER       :: fitpoints
   REAL(KIND=dp)                          :: rheavies_strength,&
                                             length, eta
   REAL(KIND=dp), DIMENSION(3)            :: box_hi, box_low
   REAL(KIND=dp), DIMENSION(:), POINTER   :: rmin_kind,&
                                             rmax_kind
   REAL(KIND=dp), DIMENSION(:), POINTER   :: range_surf 
   REAL(KIND=dp), DIMENSION(:), POINTER   :: rhs
   REAL(KIND=dp), DIMENSION(:, :),POINTER :: matrix
  END TYPE resp_type

  TYPE resp_p_type
   TYPE (resp_type), POINTER              ::  p_resp 
  END TYPE resp_p_type

CONTAINS

! *****************************************************************************
!> \brief performs resp fit and generates RESP charges
!> \param qs_env the qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE resp_fit(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'resp_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, my_per, natom, &
                                                nvar, output_unit, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    LOGICAL                                  :: failure, has_resp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(section_vals_type), POINTER         :: cons_section, input, &
                                                poisson_section, &
                                                resp_section, rest_section

    CALL timeset(routineN,handle)

    NULLIFY(logger,atomic_kind_set,cell,subsys,particles,particle_set,input,&
            resp_section,cons_section,rest_section,poisson_section,resp_env,rep_sys)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, input=input,&
            subsys=subsys,particle_set=particle_set, cell=cell, error=error)
       resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                  error=error)
       CALL section_vals_get(resp_section, explicit=has_resp, error=error)
    END IF

    IF (.NOT. failure .AND. has_resp) THEN
       logger => cp_error_get_logger(error)
       poisson_section => section_vals_get_subs_vals(input,"DFT%POISSON",error=error)
       CALL section_vals_val_get(poisson_section,"PERIODIC",i_val=my_per,error=error)
       CALL create_resp_type(resp_env, rep_sys, error) 
       !initialize the RESP fitting, get all the keywords
       CALL init_resp(qs_env,resp_env,rep_sys,input,subsys,particle_set,atomic_kind_set,&
            cell,resp_section,cons_section,rest_section,error)
       
       !print info
       CALL print_resp_parameter_info(qs_env,resp_env,rep_sys,my_per,error)
       
       CALL cp_subsys_get(subsys,particles=particles,error=error)
       natom=particles%n_els
       nvar=natom+resp_env%ncons

       ALLOCATE(ipiv(nvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       IF(.NOT.ASSOCIATED(resp_env%matrix)) THEN 
        ALLOCATE(resp_env%matrix(nvar,nvar),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ENDIF
       IF(.NOT.ASSOCIATED(resp_env%rhs)) THEN
        ALLOCATE(resp_env%rhs(nvar),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ENDIF
       ipiv =0
       resp_env%matrix = 0.0_dp
       resp_env%rhs = 0.0_dp
      
       ! calculate the matrix and the vector rhs 
       SELECT CASE (my_per)
       CASE(use_perd_none)  
        CALL calc_resp_matrix_nonper(qs_env,resp_env,atomic_kind_set, particles,cell,&
                                     resp_env%matrix,resp_env%rhs,natom,error)   
       CASE(use_perd_xyz)
        CALL calc_resp_matrix_periodic(qs_env,resp_env,rep_sys,particles,cell,natom, error)
       CASE DEFAULT
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="RESP charges only implemented for nonperiodic systems"//&
             " or XYZ periodicity!", &
             error=error, error_level=cp_failure_level)
       END SELECT
        
       output_unit=cp_print_key_unit_nr(logger,resp_section,"PRINT%PROGRAM_RUN_INFO",&
                                           extension=".resp",error=error)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T3,A,T69,I12)') "Number of potential fitting "//&
                                              "points found: ",resp_env%npoints
          WRITE(output_unit,'()')
       ENDIF
     
       !adding restraints and constraints
       CALL add_restraints_and_constraints(qs_env,resp_env,rest_section,&
            subsys,natom,cons_section,particle_set,error)

       !solve system for the values of the charges and the lagrangian multipliers
       CALL DGETRF(nvar,nvar,resp_env%matrix,nvar,ipiv,info)
       CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

       CALL DGETRS('N',nvar,1,resp_env%matrix,nvar,ipiv,resp_env%rhs,nvar,info)
       CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

       CALL print_resp_charges(qs_env,resp_env,output_unit,natom,error)
       CALL print_pot_from_resp_charges(qs_env,resp_env,particles,natom,output_unit,error)

       DEALLOCATE(ipiv, stat=stat)
       CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
       CALL resp_dealloc(resp_env,rep_sys,error)
       CALL cp_print_key_finished_output(output_unit,logger,resp_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE resp_fit

! *****************************************************************************
!> \brief creates the resp_type structure
!> \param resp_env the resp environment
!> \param rep_sys structure for repeating input sections defining fit points
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE create_resp_type(resp_env, rep_sys, error)
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_resp_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(resp_env)) CALL resp_dealloc(resp_env,rep_sys,error)
    ALLOCATE(resp_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

    IF (.NOT. failure) THEN
       NULLIFY(resp_env%matrix,&
               resp_env%fitpoints,&
               resp_env%rmin_kind,&
               resp_env%rmax_kind,&
               resp_env%rhs)

    ENDIF
     
    resp_env%equal_charges=.FALSE.
    resp_env%itc=.FALSE.
    resp_env%nonperiodic_sys=.FALSE.
    resp_env%rheavies=.FALSE.
 
    resp_env%box_hi=0.0_dp
    resp_env%box_low=0.0_dp

    resp_env%ncons=0
    resp_env%ncons_sec=0
    resp_env%nres=0
    resp_env%nrest_sec=0
    resp_env%npoints=0
    resp_env%npoints_proc=0 

  END SUBROUTINE create_resp_type

! *****************************************************************************
!> \brief deallocates the resp_type structure 
!> \param resp_env the resp environment
!> \param rep_sys structure for repeating input sections defining fit points
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE resp_dealloc(resp_env,rep_sys, error)
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'resp_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat

    IF (ASSOCIATED(resp_env)) THEN
     IF (ASSOCIATED(resp_env%matrix)) THEN
      DEALLOCATE(resp_env%matrix, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDIF
     IF (ASSOCIATED(resp_env%rhs)) THEN
      DEALLOCATE(resp_env%rhs, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDIF
     IF (ASSOCIATED(resp_env%fitpoints)) THEN
      DEALLOCATE(resp_env%fitpoints, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDIF
     IF (ASSOCIATED(resp_env%rmin_kind))THEN
      DEALLOCATE(resp_env%rmin_kind, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDIF
     IF (ASSOCIATED(resp_env%rmax_kind))THEN
      DEALLOCATE(resp_env%rmax_kind, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDIF
     DEALLOCATE(resp_env, stat=stat)
     CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    ENDIF
    IF (ASSOCIATED(rep_sys)) THEN
     DO i=1, SIZE(rep_sys)
      DEALLOCATE(rep_sys(i)%p_resp%atom_surf_list, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
      DEALLOCATE(rep_sys(i)%p_resp, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDDO
     DEALLOCATE(rep_sys, stat=stat)
     CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    ENDIF

  END SUBROUTINE resp_dealloc

! *****************************************************************************
!> \brief intializes the resp fit. Getting the parameters
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param rep_sys structure for repeating input sections defining fit points
!> \param input ...
!> \param subsys ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param cell parameters related to the simulation cell
!> \param resp_section resp section
!> \param cons_section constraints section, part of resp section
!> \param rest_section restraints section, part of resp section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE init_resp(qs_env,resp_env,rep_sys,input,subsys,particle_set,&
             atomic_kind_set,cell,resp_section,cons_section,rest_section,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: resp_section, cons_section, &
                                                rest_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_resp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nrep, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list_cons, my_stride
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: nonperiodic_section, &
                                                periodic_section

    CALL timeset(routineN,handle)

    NULLIFY(atom_list_cons, my_stride, nonperiodic_section, periodic_section)
    failure=.FALSE.

    ! get the subsections
    nonperiodic_section=>section_vals_get_subs_vals(resp_section,"NONPERIODIC_SYS",&
                                                            error=error)
    periodic_section=>section_vals_get_subs_vals(resp_section,"PERIODIC_SYS",&
                                                            error=error)
    cons_section=>section_vals_get_subs_vals(resp_section,"CONSTRAINT",&
                                                             error=error)
    rest_section=>section_vals_get_subs_vals(resp_section,"RESTRAINT",&
                                                            error=error)

    ! get and set the parameters for nonperiodic (non-surface) systems
    CALL get_parameter_nonperiodic_sys(resp_env, nonperiodic_section, cell, &
         atomic_kind_set, error)

    ! get the parameter for periodic/surface systems
    CALL section_vals_get(periodic_section, explicit=explicit, n_repetition=nrep, error=error)
    IF(explicit) THEN
     ALLOCATE(rep_sys(nrep), stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
     DO i=1, nrep
      ALLOCATE(rep_sys(i)%p_resp, stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
      NULLIFY(rep_sys(i)%p_resp%range_surf, rep_sys(i)%p_resp%atom_surf_list)
      CALL section_vals_val_get(periodic_section,"RANGE",r_vals=rep_sys(i)%p_resp%range_surf,&
                                i_rep_section=i, error=error)
      CALL section_vals_val_get(periodic_section,"LENGTH", r_val=rep_sys(i)%p_resp%length,&
                                i_rep_section=i, error=error)
      CALL section_vals_val_get(periodic_section,"SURF_DIRECTION",&
                                i_rep_section=i, i_val=rep_sys(i)%p_resp%my_fit,error=error)
      IF(ANY(rep_sys(i)%p_resp%range_surf<0.0_dp)) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Numbers in RANGE in PERIODIC_SYS cannot be negative.")
      ENDIF
      IF(rep_sys(i)%p_resp%length<=EPSILON(0.0_dp)) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Parameter LENGTH in PERIODIC_SYS has to be larger than zero.")
      ENDIF
      !list of atoms specifing the surface
      CALL build_atom_list(periodic_section,subsys,rep_sys(i)%p_resp%atom_surf_list,rep=i,error=error)
     ENDDO
    ENDIF
     

    ! get the parameters for the constraint and restraint sections    
    CALL section_vals_get(cons_section, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_get(cons_section,n_repetition=resp_env%ncons_sec,error=error)
       DO i=1,resp_env%ncons_sec
        CALL section_vals_val_get(cons_section,"EQUAL_CHARGES",l_val=resp_env%equal_charges,& 
                                  explicit=explicit,error=error)
        IF(.NOT.explicit) CYCLE
        CALL build_atom_list(cons_section,subsys,atom_list_cons,i,error=error)
        !instead of using EQUAL_CHARGES the constraint sections could be repeated
        resp_env%ncons=resp_env%ncons+SIZE(atom_list_cons)-2
        DEALLOCATE(atom_list_cons,stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
       ENDDO
    ENDIF
    CALL section_vals_get(rest_section, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_get(rest_section,n_repetition=resp_env%nrest_sec,error=error)
    ENDIF
    resp_env%ncons=resp_env%ncons+resp_env%ncons_sec
    resp_env%nres=resp_env%nres+resp_env%nrest_sec
    
    ! get the general keywords
    CALL section_vals_val_get(resp_section,"INTEGER_TOTAL_CHARGE",&
                                     l_val=resp_env%itc,error=error)
    IF (resp_env%itc) resp_env%ncons=resp_env%ncons+1

    CALL section_vals_val_get(resp_section,"RESTRAIN_HEAVIES_TO_ZERO",&
                                    l_val=resp_env%rheavies,error=error)
    IF (resp_env%rheavies) THEN
        CALL section_vals_val_get(resp_section,"RESTRAIN_HEAVIES_STRENGTH",&
                               r_val=resp_env%rheavies_strength,error=error)
    ENDIF
    CALL section_vals_val_get(resp_section,"STRIDE",i_vals=my_stride,error=error)
    CALL cp_assert(SIZE(my_stride)==1.OR.SIZE(my_stride)==3,cp_fatal_level,cp_assertion_failed,routineP,&
         "STRIDE keyword can accept only 1 (the same for X,Y,Z) or 3 values. Correct your input file."//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    IF (SIZE(my_stride)==1) THEN
       DO i = 1,3
          resp_env%stride(i) = my_stride(1)
       END DO
    ELSE
       resp_env%stride = my_stride(1:3)
    END IF
    CALL section_vals_val_get(resp_section,"WIDTH", r_val=resp_env%eta, error=error)

    CALL timestop(handle)

  END SUBROUTINE init_resp

! *****************************************************************************
!> \brief getting the parameters for nonperiodic/non-surface systems
!> \param resp_env the resp environment
!> \param nonperiodic_section input section setting parameters for sampling
!>        fitting in spheres around the atom
!> \param cell parameters related to the simulation cell
!> \param atomic_kind_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE get_parameter_nonperiodic_sys(resp_env,nonperiodic_section,cell,&
             atomic_kind_set,error)

    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(section_vals_type), POINTER         :: nonperiodic_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'get_parameter_nonperiodic_sys', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: symbol
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: ikind, j, kind_number, nkind, &
                                                nrep_rmax, nrep_rmin, stat
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: rmax, rmin
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure=.FALSE.
    nrep_rmin=0
    nrep_rmax=0
    nkind=SIZE(atomic_kind_set)

    CALL section_vals_get(nonperiodic_section, explicit=explicit, error=error)
    IF(explicit) THEN
      resp_env%nonperiodic_sys=.TRUE.
      CALL section_vals_val_get(nonperiodic_section,"RMIN",r_val=rmin,&
                                                                    error=error)
      CALL section_vals_val_get(nonperiodic_section,"RMAX",r_val=rmax,&
                                                                    error=error)
      CALL section_vals_val_get(nonperiodic_section,"RMIN_KIND", n_rep_val=nrep_rmin, error=error)
      CALL section_vals_val_get(nonperiodic_section,"RMAX_KIND", n_rep_val=nrep_rmax, error=error)
      ALLOCATE(resp_env%rmin_kind(nkind),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
      ALLOCATE(resp_env%rmax_kind(nkind),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
      resp_env%rmin_kind(:)=rmin
      resp_env%rmax_kind(:)=rmax
      DO j=1, nrep_rmin
       CALL section_vals_val_get(nonperiodic_section,"RMIN_KIND", i_rep_val=j,&
                       c_vals=tmpstringlist, error=error)
       DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind,element_symbol=symbol,kind_number=kind_number)
        IF(TRIM(tmpstringlist(2))==TRIM(symbol)) THEN
          READ(tmpstringlist(1),*) resp_env%rmin_kind(kind_number)
          resp_env%rmin_kind(kind_number)=cp_unit_to_cp2k(resp_env%rmin_kind(kind_number),&
                                         "angstrom",error=error)
        ENDIF
       ENDDO
      ENDDO
      DO j=1, nrep_rmax
       CALL section_vals_val_get(nonperiodic_section,"RMAX_KIND", i_rep_val=j,&
                       c_vals=tmpstringlist, error=error)
       DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind,element_symbol=symbol,kind_number=kind_number)
        IF(TRIM(tmpstringlist(2))==TRIM(symbol)) THEN
          READ(tmpstringlist(1),*) resp_env%rmax_kind(kind_number)
          resp_env%rmax_kind(kind_number)=cp_unit_to_cp2k(resp_env%rmax_kind(kind_number),&
                                         "angstrom",error=error)
        ENDIF
       ENDDO
      ENDDO
      
      resp_env%box_hi=(/cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3)/)
      resp_env%box_low=0.0_dp
      CALL section_vals_val_get(nonperiodic_section,"X_HI",explicit=explicit,error=error)
      IF (explicit) CALL section_vals_val_get(nonperiodic_section,"X_HI",&
                                      r_val=resp_env%box_hi(1),error=error)
      CALL section_vals_val_get(nonperiodic_section,"X_LOW",explicit=explicit,error=error)
      IF (explicit) CALL section_vals_val_get(nonperiodic_section,"X_LOW",&
                                      r_val=resp_env%box_low(1),error=error)
      CALL section_vals_val_get(nonperiodic_section,"Y_HI",explicit=explicit,error=error)
      IF (explicit) CALL section_vals_val_get(nonperiodic_section,"Y_HI",&
                                      r_val=resp_env%box_hi(2),error=error)
      CALL section_vals_val_get(nonperiodic_section,"Y_LOW",explicit=explicit,error=error)
      IF (explicit) CALL section_vals_val_get(nonperiodic_section,"Y_LOW",&
                                      r_val=resp_env%box_low(2),error=error)
      CALL section_vals_val_get(nonperiodic_section,"Z_HI",explicit=explicit,error=error)
      IF (explicit) CALL section_vals_val_get(nonperiodic_section,"Z_HI",&
                                      r_val=resp_env%box_hi(3),error=error)
      CALL section_vals_val_get(nonperiodic_section,"Z_LOW",explicit=explicit,error=error)
      IF (explicit) CALL section_vals_val_get(nonperiodic_section,"Z_LOW",&
                                      r_val=resp_env%box_low(3),error=error)
    ENDIF

  END SUBROUTINE get_parameter_nonperiodic_sys

! *****************************************************************************
!> \brief building atom lists for different sections of RESP
!> \param section input section
!> \param subsys ...
!> \param atom_list list of atoms for restraints, constraints and fit point 
!>        sampling for slab-like systems
!> \param rep input section can be repeated, this param defines for which 
!>        repetition of the input section the atom_list is built
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE build_atom_list(section,subsys,atom_list,rep,error)

    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_subsys_type), POINTER            :: subsys
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    INTEGER, INTENT(IN), OPTIONAL            :: rep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_atom_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, i, &
                                                irep, j, max_index, n_var, &
                                                num_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: indexes
    LOGICAL                                  :: failure, index_in_range

    CALL timeset(routineN,handle)
 
    NULLIFY(indexes)
    failure=.FALSE.
    irep=1
    IF(PRESENT(rep)) irep=rep 

    CALL section_vals_val_get(section,"ATOM_LIST",i_rep_section=irep,&
                                          n_rep_val=n_var,error=error)
    num_atom=0
    DO i=1,n_var
    CALL section_vals_val_get(section,"ATOM_LIST",i_rep_section=irep,&
                               i_rep_val=i,i_vals=indexes,error=error)
    num_atom=num_atom + SIZE(indexes)
    ENDDO 
    ALLOCATE(atom_list(num_atom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    atom_list=0
    num_atom=1
    DO i=1,n_var
     CALL section_vals_val_get(section,"ATOM_LIST",i_rep_section=irep,&
                                i_rep_val=i,i_vals=indexes,error=error)
     atom_list(num_atom:num_atom+SIZE(indexes)-1)=indexes(:)
     num_atom = num_atom + SIZE(indexes)
    ENDDO
    !check atom list
    num_atom=num_atom-1
    max_index=SIZE(subsys%particles%els)
    CPPrecondition(SIZE(atom_list) /= 0,cp_failure_level,routineP,error,failure)
    index_in_range=(MAXVAL(atom_list)<= max_index)&
                         .AND.(MINVAL(atom_list) > 0)
    CPPostcondition(index_in_range,cp_failure_level,routineP,error,failure)
    DO i=1,num_atom
     DO j=i+1,num_atom
      atom_a=atom_list(i)
      atom_b=atom_list(j)
      IF(atom_a==atom_b) &
      CALL stop_program(routineN,moduleN,__LINE__,&
           "There are atoms doubled in atom list for RESP.")
     ENDDO
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE build_atom_list
 
! *****************************************************************************
!> \brief build matrix and vector for nonperiodic RESP fitting
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param atomic_kind_set ...
!> \param particles ...
!> \param cell parameters related to the simulation cell
!> \param matrix coefficient matrix of the linear system of equations
!> \param rhs vector of the linear system of equations
!> \param natom number of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calc_resp_matrix_nonper(qs_env,resp_env,atomic_kind_set,particles,&
                                     cell,matrix,rhs,natom,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rhs
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_resp_matrix_nonper', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bo(2,3), gbo(2,3), handle, i, &
                                                ikind, jx, jy, jz, k, &
                                                kind_number, l, m, nkind, &
                                                now, np(3), p, stat
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: not_in_range
    REAL(KIND=dp)                            :: dh(3,3), dvol, r(3), rmax, &
                                                rmin, vec(3), vec_pbc(3), vj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dist
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_type), POINTER                   :: v_hartree_pw
    TYPE(section_vals_type), POINTER         :: input, print_key, resp_section

    CALL timeset(routineN,handle)

    NULLIFY(logger,input,particle_set,print_key,resp_section,v_hartree_pw)
    failure=.FALSE.

    IF (.NOT.cell%orthorhombic) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Nonperiodic solution for RESP charges only"//&
            " implemented for orthorhombic cells!", &
            error=error, error_level=cp_failure_level)
    END IF
    IF(.NOT.resp_env%nonperiodic_sys) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Nonperiodic solution for RESP charges (i.e. nonperiodic"//&
             " Poisson solver) can only be used with section NONPERIODC_SYS")
    ENDIF
    CALL get_qs_env(qs_env=qs_env,input=input,particle_set=particle_set,error=error)
    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                error=error)
    print_key => section_vals_get_subs_vals(resp_section,&
                                            "PRINT%V_RESP_CUBE",error=error)
    logger => cp_error_get_logger(error)

    v_hartree_pw => qs_env%ks_env%v_hartree_rspace%pw
    bo=v_hartree_pw%pw_grid%bounds_local
    gbo=v_hartree_pw%pw_grid%bounds
    np=v_hartree_pw%pw_grid%npts
    dh=v_hartree_pw%pw_grid%dh
    dvol=v_hartree_pw%pw_grid%dvol
    nkind=SIZE(atomic_kind_set)

    ALLOCATE(dist(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    ALLOCATE(not_in_range(natom,2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

    !only store fitting points when printing RESP potential
    IF (BTEST(cp_print_key_should_output(logger%iter_info,resp_section,&
        "PRINT%V_RESP_CUBE",error=error),cp_p_file)) THEN
     IF(.NOT.ASSOCIATED(resp_env%fitpoints)) THEN
        now = 1000
        ALLOCATE(resp_env%fitpoints(3,now),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
     ELSE
        now = SIZE(resp_env%fitpoints,2)
     END IF
    ENDIF

    DO jz=bo(1,3),bo(2,3)
    DO jy=bo(1,2),bo(2,2)
    DO jx=bo(1,1),bo(2,1)
       IF (.NOT.(MODULO(jz,resp_env%stride(3))==0)) CYCLE
       IF (.NOT.(MODULO(jy,resp_env%stride(2))==0)) CYCLE
       IF (.NOT.(MODULO(jx,resp_env%stride(1))==0)) CYCLE
       !bounds bo reach from -np/2 to np/2. shift of np/2 so that r(1,1,1)=(0,0,0)
       l=jx - gbo(1,1)
       k=jy - gbo(1,2)
       p=jz - gbo(1,3)
       r(3)=p*dh(3,3)+k*dh(3,2)+l*dh(3,1)
       r(2)=p*dh(2,3)+k*dh(2,2)+l*dh(2,1)
       r(1)=p*dh(1,3)+k*dh(1,2)+l*dh(1,1)
       IF (r(3)<resp_env%box_low(3).OR.r(3)>resp_env%box_hi(3)) CYCLE
       IF (r(2)<resp_env%box_low(2).OR.r(2)>resp_env%box_hi(2)) CYCLE
       IF (r(1)<resp_env%box_low(1).OR.r(1)>resp_env%box_hi(1)) CYCLE
       ! compute distance from the grid point to all atoms
       not_in_range=.FALSE.
       DO i=1,natom
          vec=r-particles%els(i)%r
          vec_pbc(1) = vec(1) - cell%hmat(1,1)*ANINT(cell%h_inv(1,1)*vec(1))
          vec_pbc(2) = vec(2) - cell%hmat(2,2)*ANINT(cell%h_inv(2,2)*vec(2))
          vec_pbc(3) = vec(3) - cell%hmat(3,3)*ANINT(cell%h_inv(3,3)*vec(3))
          dist(i)=SQRT(SUM(vec_pbc**2))
          CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind,&
               kind_number=kind_number)
          DO ikind=1,nkind
            IF(ikind==kind_number) THEN
             rmin=resp_env%rmin_kind(ikind)
             rmax=resp_env%rmax_kind(ikind)
             EXIT
            ENDIF
          ENDDO
          IF(dist(i)<rmin+1.0E-13_dp) not_in_range(i,1)=.TRUE.
          IF(dist(i)>rmax-1.0E-13_dp) not_in_range(i,2)=.TRUE.
       ENDDO
       ! check if the point is sufficiently close and  far. if OK, we can use
       ! the point for fitting, add/subtract 1.0E-13 to get rid of rounding errors when shifting atoms 
       IF(ANY(not_in_range(:,1)).OR.ALL(not_in_range(:,2))) CYCLE
       resp_env%npoints = resp_env%npoints + 1
       IF (BTEST(cp_print_key_should_output(logger%iter_info,resp_section,&
           "PRINT%V_RESP_CUBE",error=error),cp_p_file)) THEN
        IF(resp_env%npoints > now) THEN
           now = 2*now
           CALL reallocate(resp_env%fitpoints,1,3,1,now)
        ENDIF
        resp_env%fitpoints(1,resp_env%npoints) = jx
        resp_env%fitpoints(2,resp_env%npoints) = jy
        resp_env%fitpoints(3,resp_env%npoints) = jz
       ENDIF
       ! correct for the fact that v_hartree is scaled by dvol, and has the opposite sign
       IF (qs_env%qmmm) THEN
          ! If it's a QM/MM run let's remove the contribution of the MM potential out of the Hartree pot
          vj=-v_hartree_pw%cr3d(jx,jy,jz)/dvol+qs_env%ks_qmmm_env%v_qmmm_rspace%pw%cr3d(jx,jy,jz)
       ELSE                                                                                       
          vj=-v_hartree_pw%cr3d(jx,jy,jz)/dvol
       END IF
       dist=1.0_dp/dist

       DO i=1,natom
        DO m=1,natom
           matrix(m,i)=matrix(m,i)+2.0_dp*dist(i)*dist(m)
        ENDDO
        rhs(i)=rhs(i)+2.0_dp*vj*dist(i)
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    resp_env%npoints_proc=resp_env%npoints
    CALL mp_sum(resp_env%npoints,v_hartree_pw%pw_grid%para%group)
    CALL mp_sum(matrix,v_hartree_pw%pw_grid%para%group)
    CALL mp_sum(rhs,v_hartree_pw%pw_grid%para%group)
    !weighted sum
    matrix=matrix/resp_env%npoints
    rhs=rhs/resp_env%npoints

    DEALLOCATE(dist,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    DEALLOCATE(not_in_range,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)

    CALL timestop(handle)

  END SUBROUTINE calc_resp_matrix_nonper
 
! *****************************************************************************
!> \brief build matrix and vector for periodic RESP fitting
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param rep_sys structure for repeating input sections defining fit points
!> \param particles ...
!> \param cell parameters related to the simulation cell
!> \param natom number of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calc_resp_matrix_periodic(qs_env,resp_env,rep_sys,particles,cell,&
                                       natom,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_resp_matrix_periodic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ip, j, jx, jy, jz, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: normalize_factor
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vpot
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_ga, va_gspace, va_rspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(pw_env, para_env,auxbas_pw_pool,poisson_env)

    IF(.NOT.ALL(cell%perd/=0)) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Periodic solution for RESP (with periodic Poisson solver)"//&
             " can only be obtained with a cell that has XYZ periodicity")
    ENDIF

    CALL get_qs_env(qs_env, pw_env=pw_env,para_env=para_env,&
                      error=error)   

    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                      poisson_env=poisson_env, error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                             rho_ga%pw,&
                             use_data=COMPLEXDATA1D,&
                             in_space=RECIPROCALSPACE,&
                             error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                             va_gspace%pw,&
                             use_data=COMPLEXDATA1D,&
                             in_space=RECIPROCALSPACE,&
                             error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                             va_rspace%pw,&
                             use_data=REALDATA3D,&
                             in_space=REALSPACE,&
                             error=error)

    !get fitting points and store them in resp_env%fitpoints
    CALL get_fitting_points(qs_env,resp_env,rep_sys,particles=particles,&
         cell=cell,error=error)
    ALLOCATE(vpot(resp_env%npoints,natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    normalize_factor=SQRT((resp_env%eta/pi)**3)

    DO i=1,natom
     !collocate gaussian for each atom
      CALL pw_zero(rho_ga%pw, error=error)
      CALL calculate_rho_resp_single(rho_ga,qs_env,resp_env%eta,i,error)
     !calculate potential va and store the part needed for fitting in vpot
      CALL pw_zero(va_gspace%pw, error=error)
      CALL pw_poisson_solve(poisson_env,rho_ga%pw,vhartree=va_gspace%pw,error=error)
      CALL pw_zero(va_rspace%pw, error=error)
      CALL pw_transfer(va_gspace%pw,va_rspace%pw,error=error)
      CALL pw_scale(va_rspace%pw,normalize_factor,error=error)
      DO ip=1,resp_env%npoints
         jx = resp_env%fitpoints(1,ip)
         jy = resp_env%fitpoints(2,ip)
         jz = resp_env%fitpoints(3,ip)
         vpot(ip,i) = va_rspace%pw%cr3d(jx,jy,jz)
      END DO
    ENDDO

    CALL pw_release(va_gspace%pw,error=error)
    CALL pw_release(va_rspace%pw,error=error)
    CALL pw_release(rho_ga%pw,error=error)

    DO i=1,natom
      DO j=1,natom
      ! calculate matrix
         resp_env%matrix(i,j)=resp_env%matrix(i,j) + 2.0_dp*SUM(vpot(:,i)*vpot(:,j))
      ENDDO
      ! calculate vector resp_env%rhs
      CALL calculate_rhs(qs_env,resp_env,resp_env%rhs(i),vpot(:,i),error=error)
    ENDDO

    DEALLOCATE(vpot,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    
    resp_env%npoints_proc=resp_env%npoints
    CALL mp_sum(resp_env%npoints,para_env%group)
    CALL mp_sum(resp_env%matrix,para_env%group)
    CALL mp_sum(resp_env%rhs,para_env%group)
    !weighted sum
    resp_env%matrix=resp_env%matrix/resp_env%npoints
    resp_env%rhs=resp_env%rhs/resp_env%npoints

    CALL timestop(handle)

  END SUBROUTINE calc_resp_matrix_periodic

! *****************************************************************************
!> \brief get RESP fitting points for the periodic fitting
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param rep_sys structure for repeating input sections defining fit points
!> \param particles ...
!> \param cell parameters related to the simulation cell
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE get_fitting_points(qs_env,resp_env,rep_sys,particles,cell,error)
 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_fitting_points', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2,3), gbo(2,3), handle, i, iatom, ikind, in_x, in_y, in_z, &
      jx, jy, jz, k, kind_number, l, m, natom, nkind, now, output_unit, p, &
      stat
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: not_in_range
    REAL(KIND=dp)                            :: dh(3,3), r(3), rmax, rmin, &
                                                vec_pbc(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dist
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_type), POINTER                   :: v_hartree_pw
    TYPE(section_vals_type), POINTER         :: input, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(atomic_kind_set,v_hartree_pw,para_env,logger,input,particle_set,&
            resp_section)   
 
    CALL get_qs_env(qs_env,input=input,particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,para_env=para_env,&
                    error=error)

    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                error=error)
    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,resp_section,&
                                    "PRINT%COORD_FIT_POINTS",&
                                     extension=".xyz",&
                                     file_status="REPLACE",&
                                     file_action="WRITE",&
                                     file_form="FORMATTED",& 
                                     error=error)

    v_hartree_pw => qs_env%ks_env%v_hartree_rspace%pw
    bo=v_hartree_pw%pw_grid%bounds_local
    gbo=v_hartree_pw%pw_grid%bounds
    dh=v_hartree_pw%pw_grid%dh
    natom=SIZE(particles%els)
    nkind=SIZE(atomic_kind_set)

    IF(.NOT.ASSOCIATED(resp_env%fitpoints)) THEN
       now = 1000
       ALLOCATE(resp_env%fitpoints(3,now),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    ELSE
       now = SIZE(resp_env%fitpoints,2)
    END IF

    ALLOCATE(dist(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    ALLOCATE(not_in_range(natom,2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

    !every proc gets another bo, grid is distributed
    DO jz=bo(1,3),bo(2,3)
     IF (.NOT.(MODULO(jz,resp_env%stride(3))==0)) CYCLE
     DO jy=bo(1,2),bo(2,2)
      IF (.NOT.(MODULO(jy,resp_env%stride(2))==0)) CYCLE
      DO jx=bo(1,1),bo(2,1)
       IF (.NOT.(MODULO(jx,resp_env%stride(1))==0)) CYCLE
       !bounds gbo reach from -np/2 to np/2. shift of np/2 so that r(1,1,1)=(0,0,0)
       l=jx - gbo(1,1)
       k=jy - gbo(1,2)
       p=jz - gbo(1,3)
       r(3)=p*dh(3,3)+k*dh(3,2)+l*dh(3,1)
       r(2)=p*dh(2,3)+k*dh(2,2)+l*dh(2,1)
       r(1)=p*dh(1,3)+k*dh(1,2)+l*dh(1,1)
       IF(resp_env%nonperiodic_sys) THEN
         not_in_range=.FALSE.
         DO m=1,natom
            vec_pbc = pbc(r,particles%els(m)%r,cell)
            dist(m)=SQRT(SUM(vec_pbc**2))
            CALL get_atomic_kind(atomic_kind=particle_set(m)%atomic_kind,&
                 kind_number=kind_number)
            DO ikind=1,nkind
              IF(ikind==kind_number) THEN
               rmin=resp_env%rmin_kind(ikind)
               rmax=resp_env%rmax_kind(ikind)
               EXIT
              ENDIF
            ENDDO
            IF(dist(m)<rmin+1.0E-13_dp) not_in_range(m,1)=.TRUE.
            IF(dist(m)>rmax-1.0E-13_dp) not_in_range(m,2)=.TRUE.
         ENDDO
         IF(ANY(not_in_range(:,1)).OR.ALL(not_in_range(:,2))) CYCLE
       ELSE
         DO i=1,SIZE(rep_sys)
          DO m=1,SIZE(rep_sys(i)%p_resp%atom_surf_list)
             in_z=0
             in_y=0
             in_x=0
             iatom=rep_sys(i)%p_resp%atom_surf_list(m)
             SELECT CASE(rep_sys(i)%p_resp%my_fit)
             CASE(do_resp_x_dir, do_resp_y_dir, do_resp_z_dir)
              vec_pbc = pbc(particles%els(iatom)%r,r,cell)
             CASE(do_resp_minus_x_dir, do_resp_minus_y_dir, do_resp_minus_z_dir)
              vec_pbc = pbc(r,particles%els(iatom)%r,cell)
             END SELECT
             SELECT CASE(rep_sys(i)%p_resp%my_fit)
             !subtract 1.0E-13 to get rid of rounding errors when shifting atoms 
             CASE(do_resp_x_dir, do_resp_minus_x_dir)
              IF(ABS(vec_pbc(3))<rep_sys(i)%p_resp%length-1.0E-13_dp)       in_z=1
              IF(ABS(vec_pbc(2))<rep_sys(i)%p_resp%length-1.0E-13_dp)       in_y=1
              IF(vec_pbc(1)>rep_sys(i)%p_resp%range_surf(1)+1.0E-13_dp.AND.&              
                  vec_pbc(1)<rep_sys(i)%p_resp%range_surf(2)-1.0E-13_dp)    in_x=1 
             CASE(do_resp_y_dir,do_resp_minus_y_dir)
              IF(ABS(vec_pbc(3))<rep_sys(i)%p_resp%length-1.0E-13_dp)       in_z=1
              IF(vec_pbc(2)>rep_sys(i)%p_resp%range_surf(1)+1.0E-13_dp.AND.&              
                  vec_pbc(2)<rep_sys(i)%p_resp%range_surf(2)-1.0E-13_dp)    in_y=1 
              IF(ABS(vec_pbc(1))<rep_sys(i)%p_resp%length-1.0E-13_dp)       in_x=1
             CASE(do_resp_z_dir, do_resp_minus_z_dir)
              IF(vec_pbc(3)>rep_sys(i)%p_resp%range_surf(1)+1.0E-13_dp.AND.&              
                  vec_pbc(3)<rep_sys(i)%p_resp%range_surf(2)-1.0E-13_dp)    in_z=1 
              IF(ABS(vec_pbc(2))<rep_sys(i)%p_resp%length-1.0E-13_dp)       in_y=1
              IF(ABS(vec_pbc(1))<rep_sys(i)%p_resp%length-1.0E-13_dp)       in_x=1
             END SELECT
             IF(in_z*in_y*in_x==1) EXIT
          ENDDO
          IF(in_z*in_y*in_x==1) EXIT
         ENDDO
         IF(in_z*in_y*in_x==0) CYCLE
       ENDIF
       resp_env%npoints=resp_env%npoints+1
       IF(resp_env%npoints > now) THEN
          now = 2*now
          CALL reallocate(resp_env%fitpoints,1,3,1,now)
       ENDIF
       resp_env%fitpoints(1,resp_env%npoints) = jx
       resp_env%fitpoints(2,resp_env%npoints) = jy
       resp_env%fitpoints(3,resp_env%npoints) = jz
      ENDDO
     ENDDO
    ENDDO

    !print fitting points to file if requested 
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
            resp_section,"PRINT%COORD_FIT_POINTS",error=error),&
            cp_p_file))THEN
     CALL print_fitting_points(qs_env,resp_env,dh,output_unit,gbo,error)
    ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,resp_section,&
                       "PRINT%COORD_FIT_POINTS", error=error)
    DEALLOCATE(dist,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    DEALLOCATE(not_in_range,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)

    CALL timestop(handle)
 
  END SUBROUTINE get_fitting_points

! *****************************************************************************
!> \brief calculate vector rhs
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param rhs vector 
!> \param vpot single gaussian potential
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calculate_rhs(qs_env,resp_env,rhs,vpot,error)
 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    REAL(KIND=dp), INTENT(INOUT)             :: rhs
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vpot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rhs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ip, jx, jy, jz
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dvol, vj
    TYPE(pw_type), POINTER                   :: v_hartree_pw

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(v_hartree_pw)
    v_hartree_pw => qs_env%ks_env%v_hartree_rspace%pw
    dvol=v_hartree_pw%pw_grid%dvol
    !multiply v_hartree and va_rspace and calculate the vector rhs
    !taking into account that v_hartree has opposite site; remove v_qmmm       
    DO ip=1,resp_env%npoints
       jx = resp_env%fitpoints(1,ip)
       jy = resp_env%fitpoints(2,ip)
       jz = resp_env%fitpoints(3,ip)
       vj = -v_hartree_pw%cr3d(jx,jy,jz)/dvol
       IF (qs_env%qmmm) THEN
         !taking into account that v_qmmm has also opposite sign
         vj = vj + qs_env%ks_qmmm_env%v_qmmm_rspace%pw%cr3d(jx,jy,jz) 
       ENDIF
       rhs = rhs + 2.0_dp*vj*vpot(ip)
    ENDDO
 
    CALL timestop(handle)
 
  END SUBROUTINE calculate_rhs

! *****************************************************************************
!> \brief print the atom coordinates and the coordinates of the fitting points 
!>        to an xyz file
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param dh incremental cell matrix
!> \param output_unit ...
!> \param gbo lower and upper bounds for the grid points
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_fitting_points(qs_env,resp_env,dh,output_unit,gbo,error)
  
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    REAL(KIND=dp), INTENT(IN)                :: dh(3,3)
    INTEGER, INTENT(IN)                      :: output_unit, gbo(2,3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_fitting_points', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, i, iatom, ip, jx, jy, &
                                                jz, k, l, my_pos, p, req(6), &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: tmp_npoints, tmp_size
    INTEGER, DIMENSION(:, :), POINTER        :: tmp_points
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: conv, r(3)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: input, print_key, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(para_env,input,logger,resp_section,print_key,particle_set,tmp_size,&
            tmp_points,tmp_npoints)   
 
    CALL get_qs_env(qs_env, input=input, para_env=para_env,&
                       particle_set=particle_set,error=error)
 
    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                error=error)
    print_key => section_vals_get_subs_vals(resp_section,&
                                            "PRINT%COORD_FIT_POINTS",&
                                                error=error)
    logger => cp_error_get_logger(error)
    conv=cp_unit_from_cp2k(1.0_dp,"angstrom",error=error)

    IF(output_unit>0) THEN
     filename=cp_print_key_generate_filename(logger,&
              print_key, extension=".xyz",&
              my_local=.FALSE.,error=error)
     WRITE(unit=output_unit,FMT="(I12,A,/)") SIZE(particle_set),' + nr fit points'
     DO iatom=1,SIZE(particle_set)
      CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                           element_symbol=element_symbol)
      WRITE(UNIT=output_unit,FMT="(A,1X,3F10.5)") element_symbol,&
                                              particle_set(iatom)%r(1:3)*conv
     ENDDO
     !printing points of proc which is doing the output (should be proc 0)
     DO ip=1,resp_env%npoints
      jx=resp_env%fitpoints(1,ip)     
      jy=resp_env%fitpoints(2,ip)     
      jz=resp_env%fitpoints(3,ip)
      l=jx - gbo(1,1)
      k=jy - gbo(1,2)
      p=jz - gbo(1,3) 
      r(3)=p*dh(3,3)+k*dh(3,2)+l*dh(3,1) 
      r(2)=p*dh(2,3)+k*dh(2,2)+l*dh(2,1)
      r(1)=p*dh(1,3)+k*dh(1,2)+l*dh(1,1)  
      r(:)=r(:)*conv 
      WRITE(UNIT=output_unit,FMT="(A,2X,3F10.5)") "X", r(1), r(2), r(3)
     ENDDO 
    ENDIF

    ALLOCATE(tmp_size(1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    ALLOCATE(tmp_npoints(1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

    !sending data of all other prosc to proc which makes the output (proc 0) 
    IF(output_unit>0) THEN
      my_pos=para_env%mepos
      DO i=1,para_env%num_pe
       IF(my_pos==i-1) CYCLE
       CALL mp_irecv(msgout=tmp_size,source=i-1,comm=para_env%group,&
            request=req(1))
       CALL mp_wait(req(1))
       ALLOCATE(tmp_points(3,tmp_size(1)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       CALL mp_irecv(msgout=tmp_points,source=i-1,comm=para_env%group,&
            request=req(3))
       CALL mp_wait(req(3))
       CALL mp_irecv(msgout=tmp_npoints,source=i-1,comm=para_env%group,&
            request=req(5))
       CALL mp_wait(req(5))
       DO ip=1,tmp_npoints(1) 
        jx=tmp_points(1,ip)     
        jy=tmp_points(2,ip)     
        jz=tmp_points(3,ip)
        l=jx - gbo(1,1)
        k=jy - gbo(1,2)
        p=jz - gbo(1,3) 
        r(3)=p*dh(3,3)+k*dh(3,2)+l*dh(3,1) 
        r(2)=p*dh(2,3)+k*dh(2,2)+l*dh(2,1)
        r(1)=p*dh(1,3)+k*dh(1,2)+l*dh(1,1)  
        r(:)=r(:)*conv 
        WRITE(UNIT=output_unit,FMT="(A,2X,3F10.5)") "X", r(1), r(2),r(3)
       ENDDO 
       DEALLOCATE(tmp_points,stat=stat)
       CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
      ENDDO
    ELSE
     tmp_size(1)=SIZE(resp_env%fitpoints,2)
     !para_env%source should be 0
     CALL mp_isend(msgin=tmp_size,dest=para_env%source,comm=para_env%group,&
          request=req(2))
     CALL mp_wait(req(2))
     CALL mp_isend(msgin=resp_env%fitpoints,dest=para_env%source,comm=para_env%group,&
          request=req(4))
     CALL mp_wait(req(4))
     tmp_npoints(1)=resp_env%npoints
     CALL mp_isend(msgin=tmp_npoints,dest=para_env%source,comm=para_env%group,&
          request=req(6))
     CALL mp_wait(req(6))
    ENDIF

    DEALLOCATE(tmp_size,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    DEALLOCATE(tmp_npoints,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)

    CALL timestop(handle)
 
  END SUBROUTINE print_fitting_points

! *****************************************************************************
!> \brief add restraints and constraints
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param rest_section input section for restraints
!> \param subsys ...
!> \param natom number of atoms
!> \param cons_section input section for constraints
!> \param particle_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE add_restraints_and_constraints(qs_env,resp_env,rest_section,&
                             subsys,natom,cons_section,particle_set,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(section_vals_type), POINTER         :: rest_section
    TYPE(cp_subsys_type), POINTER            :: subsys
    INTEGER, INTENT(IN)                      :: natom
    TYPE(section_vals_type), POINTER         :: cons_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'add_restraints_and_constraints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, k, m, ncons_v, &
                                                stat, z
    INTEGER, DIMENSION(:), POINTER           :: atom_list_cons, atom_list_res
    LOGICAL                                  :: explicit_coeff, failure
    REAL(KIND=dp)                            :: my_atom_coef(2), strength, &
                                                TARGET
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atom_coef

    CALL timeset(routineN,handle)

    NULLIFY(atom_coef,atom_list_res,atom_list_cons)
    failure=.FALSE.
 
    ! add the restraints
    DO i=1,resp_env%nrest_sec
       CALL section_vals_val_get(rest_section,"TARGET",i_rep_section=i,r_val=TARGET,error=error)
       CALL section_vals_val_get(rest_section,"STRENGTH",i_rep_section=i,r_val=strength,error=error)
       CALL build_atom_list(rest_section,subsys,atom_list_res,i,error)
       CALL section_vals_val_get(rest_section,"ATOM_COEF",i_rep_section=i,explicit=explicit_coeff,&
            error=error)
       IF (explicit_coeff) THEN
        CALL section_vals_val_get(rest_section,"ATOM_COEF",i_rep_section=i,r_vals=atom_coef,&
             error=error)
        CPPrecondition(SIZE(atom_list_res)==SIZE(atom_coef),cp_failure_level,routineP,error,failure)
       ENDIF
       DO m=1,SIZE(atom_list_res)
         IF (explicit_coeff) THEN
          DO k=1,SIZE(atom_list_res)
             resp_env%matrix(atom_list_res(m),atom_list_res(k))=&
                                     resp_env%matrix(atom_list_res(m),atom_list_res(k))+ &
                                     atom_coef(m)*atom_coef(k)*2.0_dp*strength
          ENDDO
          resp_env%rhs(atom_list_res(m))=resp_env%rhs(atom_list_res(m))+&
                                       2.0_dp*TARGET*strength*atom_coef(m)
         ELSE
          resp_env%matrix(atom_list_res(m),atom_list_res(m))=&
              resp_env%matrix(atom_list_res(m),atom_list_res(m))+ &
               2.0_dp*strength
          resp_env%rhs(atom_list_res(m))=resp_env%rhs(atom_list_res(m))+&
                                       2.0_dp*TARGET*strength
         ENDIF
       ENDDO
       DEALLOCATE(atom_list_res,stat=stat)
       CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    ENDDO

    ! if heavies are restrained to zero, add these as well
    IF (resp_env%rheavies) THEN
       DO i=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind,z=z)
          IF (z.NE.1) THEN
             resp_env%matrix(i,i)=resp_env%matrix(i,i)+2.0_dp*resp_env%rheavies_strength
          ENDIF
       ENDDO
    ENDIF

    ! add the constraints
    ncons_v=0
    ncons_v=ncons_v+natom
    IF (resp_env%itc) THEN
       ncons_v=ncons_v+1
       resp_env%matrix(1:natom,ncons_v)=1.0_dp
       resp_env%matrix(ncons_v,1:natom)=1.0_dp
       resp_env%rhs(ncons_v)=qs_env%dft_control%charge
    ENDIF

    DO i=1,resp_env%ncons_sec
       CALL build_atom_list(cons_section,subsys,atom_list_cons,i,error)
       IF(.NOT.resp_env%equal_charges) THEN
         ncons_v=ncons_v+1
         CALL section_vals_val_get(cons_section,"ATOM_COEF",i_rep_section=i,r_vals=atom_coef,error=error)
         CALL section_vals_val_get(cons_section,"TARGET",i_rep_section=i,r_val=TARGET,error=error)
         CPPrecondition(SIZE(atom_list_cons)==SIZE(atom_coef),cp_failure_level,routineP,error,failure)
         DO m=1,SIZE(atom_list_cons)
            resp_env%matrix(atom_list_cons(m),ncons_v)=atom_coef(m)
            resp_env%matrix(ncons_v,atom_list_cons(m))=atom_coef(m)
         ENDDO
         resp_env%rhs(ncons_v)=TARGET
       ELSE
         my_atom_coef(1)=1.0_dp
         my_atom_coef(2)=-1.0_dp
         DO k=2,SIZE(atom_list_cons)
            ncons_v=ncons_v+1
            resp_env%matrix(atom_list_cons(1),ncons_v)=my_atom_coef(1)
            resp_env%matrix(ncons_v,atom_list_cons(1))=my_atom_coef(1)
            resp_env%matrix(atom_list_cons(k),ncons_v)=my_atom_coef(2)
            resp_env%matrix(ncons_v,atom_list_cons(k))=my_atom_coef(2)
            resp_env%rhs(ncons_v)=0.0_dp
         ENDDO
       ENDIF
       DEALLOCATE(atom_list_cons,stat=stat)
       CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    ENDDO
    CALL timestop(handle)

  END SUBROUTINE add_restraints_and_constraints

! *****************************************************************************
!> \brief print input information
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param rep_sys structure for repeating input sections defining fit points 
!> \param my_per ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_resp_parameter_info(qs_env,resp_env,rep_sys,my_per,error)
   
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(resp_p_type), DIMENSION(:), POINTER :: rep_sys
    INTEGER, INTENT(IN)                      :: my_per
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_resp_parameter_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, output_unit
    REAL(KIND=dp)                            :: conv, eta_conv
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: input, resp_section

    CALL timeset(routineN,handle)
    NULLIFY(logger,input,resp_section)

    CALL get_qs_env(qs_env, input=input, error=error)
    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                               error=error)
    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,resp_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".resp",error=error)

    conv=cp_unit_from_cp2k(1.0_dp,"angstrom",error=error)
    IF(.NOT.my_per==use_perd_none) THEN
     eta_conv=cp_unit_from_cp2k(resp_env%eta,"angstrom",power=-2,error=error)
    ENDIF

    IF (output_unit>0) THEN
     WRITE(output_unit,'(/,1X,A,/)')      "STARTING RESP FIT"
     IF(.NOT.resp_env%equal_charges) THEN
      WRITE(output_unit,'(T3,A,T75,I6)') "Number of explicit constraints: ",resp_env%ncons_sec
     ELSE
      IF(resp_env%itc) THEN
       WRITE(output_unit,'(T3,A,T75,I6)') "Number of explicit constraints: ",resp_env%ncons-1
      ELSE
       WRITE(output_unit,'(T3,A,T75,I6)') "Number of explicit constraints: ",resp_env%ncons
      ENDIF  
     ENDIF
     WRITE(output_unit,'(T3,A,T75,I6)') "Number of explicit restraints: ",resp_env%nrest_sec
     WRITE(output_unit,'(T3,A,T80,A)')  "Constrain total charge ",MERGE("T","F",resp_env%itc)
     WRITE(output_unit,'(T3,A,T80,A)')  "Restrain heavy atoms ",MERGE("T","F",resp_env%rheavies)
     IF (resp_env%rheavies) THEN
        WRITE(output_unit,'(T3,A,T71,F10.6)') "Heavy atom restraint strength: ",&
                                                                    resp_env%rheavies_strength
     ENDIF
     WRITE(output_unit,'(T3,A,T66,3I5)') "Stride: ",resp_env%stride
     IF(resp_env%nonperiodic_sys) THEN
       IF(ALL(resp_env%rmin_kind==resp_env%rmin_kind(1)))& 
        WRITE(output_unit,'(T3,A,T71,F10.5)') "Rmin [angstrom]: ",resp_env%rmin_kind(1)*conv
       IF(ALL(resp_env%rmax_kind==resp_env%rmax_kind(1)))& 
        WRITE(output_unit,'(T3,A,T71,F10.5)')  "Rmax [angstrom]: ",resp_env%rmax_kind(1)*conv
       WRITE(output_unit,'(T3,A,T51,3F10.5)') "Box min [angstrom]: ",resp_env%box_low(1:3)*conv
       WRITE(output_unit,'(T3,A,T51,3F10.5)') "Box max [angstrom]: ",resp_env%box_hi(1:3)*conv
     ELSE
       WRITE(output_unit,'(2X,A,F10.5)')  "Index of atoms defining the surface: "
       DO i=1, SIZE(rep_sys)
        IF(i>1.AND.ALL(rep_sys(i)%p_resp%atom_surf_list==rep_sys(1)%p_resp%atom_surf_list)) EXIT
        WRITE(output_unit,'(7X,10I6)')  rep_sys(i)%p_resp%atom_surf_list
       ENDDO
       DO i=1, SIZE(rep_sys)
        IF(i>1.AND.ALL(rep_sys(i)%p_resp%range_surf==rep_sys(1)%p_resp%range_surf)) EXIT
        WRITE(output_unit,'(T3,A,T61,2F10.5)')&
                           "Range for sampling above the surface [angstrom]:",&
                            rep_sys(i)%p_resp%range_surf(1:2)*conv
       ENDDO
       DO i=1, SIZE(rep_sys)
        IF(i>1.AND.rep_sys(i)%p_resp%length==rep_sys(1)%p_resp%length) EXIT
        WRITE(output_unit,'(T3,A,T71,F10.5)')  "Length of sampling box above each"//&
                           " surface atom [angstrom]: ",rep_sys(i)%p_resp%length*conv
       ENDDO
     ENDIF
     IF(.NOT.my_per==use_perd_none) THEN
       WRITE(output_unit,'(T3,A,T71,F10.5)')  "Width of Gaussian charge"//&
                                              " distribution [angstrom^-2]: ", eta_conv
     ENDIF
     CALL m_flush(output_unit)
    ENDIF
    CALL cp_print_key_finished_output(output_unit,logger,resp_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

    CALL timestop(handle)

  END SUBROUTINE print_resp_parameter_info

! *****************************************************************************
!> \brief print RESP charges to an extra file or to the normal output file
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param output_runinfo ...
!> \param natom number of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_resp_charges(qs_env,resp_env,output_runinfo,natom,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    INTEGER, INTENT(IN)                      :: output_runinfo, natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_resp_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, output_file
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: input, print_key, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(particle_set,qs_kind_set,input,logger,resp_section,print_key)
 
    CALL get_qs_env(qs_env,input=input,particle_set=particle_set,&
                    qs_kind_set=qs_kind_set, error=error)

    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                error=error)
    print_key => section_vals_get_subs_vals(resp_section,&
                                            "PRINT%RESP_CHARGES_TO_FILE",&
                                                error=error)
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              resp_section,"PRINT%RESP_CHARGES_TO_FILE",error=error),&
              cp_p_file)) THEN
     output_file=cp_print_key_unit_nr(logger,resp_section,&
                                     "PRINT%RESP_CHARGES_TO_FILE",&
                                      extension=".resp",&
                                      file_status="REPLACE",&
                                      file_action="WRITE",&
                                      file_form="FORMATTED",& 
                                      error=error)
     IF(output_file>0) THEN
      filename = cp_print_key_generate_filename(logger,&
                 print_key, extension=".resp", &
                 my_local=.FALSE.,error=error)
     CALL print_atomic_charges(particle_set,qs_kind_set,output_file,title="RESP charges:",&
                                             atomic_charges=resp_env%rhs(1:natom))
     IF(output_runinfo>0) WRITE(output_runinfo,'(2X,A,/)')  "PRINTED RESP CHARGES TO FILE"
     ENDIF
 
     CALL cp_print_key_finished_output(output_file,logger,resp_section,&
                       "PRINT%RESP_CHARGES_TO_FILE", error=error)
    ELSE
     CALL print_atomic_charges(particle_set,qs_kind_set,output_runinfo,title="RESP charges:",&
                                             atomic_charges=resp_env%rhs(1:natom))
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE print_resp_charges

! *****************************************************************************
!> \brief print potential generated by RESP charges to file
!> \param qs_env the qs environment
!> \param resp_env the resp environment
!> \param particles ...
!> \param natom number of atoms
!> \param output_runinfo ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_pot_from_resp_charges(qs_env,resp_env,particles,natom,output_runinfo,error)
 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: natom, output_runinfo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_pot_from_resp_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: my_pos_cube
    INTEGER                                  :: handle, ip, jx, jy, jz, &
                                                unit_nr
    LOGICAL                                  :: append_cube, failure
    REAL(KIND=dp)                            :: normalize_factor, rms, rrms, &
                                                sum_diff, sum_hartree, udvol
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: aux_r, rho_resp, &
                                                v_resp_gspace, v_resp_rspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: input, print_key, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(auxbas_pw_pool,logger,pw_env,poisson_env,input,print_key,&
            para_env,resp_section)
    normalize_factor=SQRT((resp_env%eta/pi)**3)
    CALL get_qs_env(qs_env=qs_env,input=input,para_env=para_env,&
                    pw_env=pw_env,error=error)
    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                error=error)
    print_key => section_vals_get_subs_vals(resp_section,&
                                            "PRINT%V_RESP_CUBE",error=error)
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,resp_section,&
         "PRINT%V_RESP_CUBE",error=error),cp_p_file)) THEN
     ! calculate potential generated from RESP charges
     CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,&
                       poisson_env=poisson_env,error=error)

     CALL pw_pool_create_pw(auxbas_pw_pool,&
                            rho_resp%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,&
                            error=error)
     CALL pw_pool_create_pw(auxbas_pw_pool,&
                            v_resp_gspace%pw, &
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,&
                            error=error)
     CALL pw_pool_create_pw(auxbas_pw_pool,&
                            v_resp_rspace%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,&
                            error=error)

     CALL pw_zero(rho_resp%pw,error=error)
     CALL calculate_rho_resp_all(rho_resp,resp_env%rhs,natom,&
                                              resp_env%eta,qs_env,error)
     CALL pw_zero(v_resp_gspace%pw, error=error)
     CALL pw_poisson_solve(poisson_env,rho_resp%pw,&
                          vhartree=v_resp_gspace%pw,error=error)
     CALL pw_zero(v_resp_rspace%pw, error=error)
     CALL pw_transfer(v_resp_gspace%pw,v_resp_rspace%pw,error=error)
     CALL pw_scale(v_resp_rspace%pw,v_resp_rspace%pw%pw_grid%dvol,&
                   error=error) 
     CALL pw_scale(v_resp_rspace%pw,-normalize_factor,&
                   error=error)  
     CALL pw_release(v_resp_gspace%pw,error=error)
     CALL pw_release(rho_resp%pw,error=error)
     !now print the v_resp_rspace%pw to a cube file
     CALL pw_pool_create_pw(auxbas_pw_pool,aux_r%pw,&
                             use_data = REALDATA3D,&
                             in_space = REALSPACE, error=error)
     append_cube=section_get_lval(resp_section,"PRINT%V_RESP_CUBE%APPEND",error=error)
     my_pos_cube="REWIND"
     IF(append_cube) THEN
       my_pos_cube="APPEND"
     END IF
     unit_nr=cp_print_key_unit_nr(logger,resp_section,&
                           "PRINT%V_RESP_CUBE",&
                           extension=".cube",&
                           file_position=my_pos_cube,error=error)
     udvol = 1.0_dp/v_resp_rspace%pw%pw_grid%dvol
     CALL pw_copy(v_resp_rspace%pw,aux_r%pw, error=error)
     CALL pw_scale(aux_r%pw,udvol,error=error)
     CALL cp_pw_to_cube(aux_r%pw,unit_nr,"RESP POTENTIAL",particles=particles,&
               stride=section_get_ivals(resp_section,&
               "PRINT%V_RESP_CUBE%STRIDE",error=error),&
               error=error)
     CALL cp_print_key_finished_output(unit_nr,logger,resp_section,&
          "PRINT%V_RESP_CUBE",error=error)
     CALL pw_pool_give_back_pw(auxbas_pw_pool,aux_r%pw, error=error)

    !RMS and RRMS
    sum_diff=0.0_dp
    sum_hartree=0.0_dp
    rms=0.0_dp
    rrms=0.0_dp
    DO ip=1,resp_env%npoints_proc 
     jx = resp_env%fitpoints(1,ip)
     jy = resp_env%fitpoints(2,ip)
     jz = resp_env%fitpoints(3,ip)
     sum_diff=sum_diff+(qs_env%ks_env%v_hartree_rspace%pw%cr3d(jx,jy,jz)-&
              v_resp_rspace%pw%cr3d(jx,jy,jz))**2
     sum_hartree=sum_hartree+qs_env%ks_env%v_hartree_rspace%pw%cr3d(jx,jy,jz)**2 
    ENDDO
     CALL mp_sum(sum_diff,para_env%group)
     CALL mp_sum(sum_hartree,para_env%group)
     rms=SQRT(sum_diff/resp_env%npoints)
     rrms=SQRT(sum_diff/sum_hartree)
     IF(output_runinfo > 0) THEN
      WRITE(output_runinfo,'(2X,A,T69,ES12.5)') "Root-mean-square (RMS) "//&
                                                "error of RESP fit:",rms  
      WRITE(output_runinfo,'(2X,A,T69,ES12.5,/)') "Relative root-mean-square "//&
                                                  "(RRMS) error of RESP fit:",rrms  
     ENDIF
     CALL pw_release(v_resp_rspace%pw,error=error)
    ENDIF
 
    CALL timestop(handle)

  END SUBROUTINE print_pot_from_resp_charges
END MODULE qs_resp
