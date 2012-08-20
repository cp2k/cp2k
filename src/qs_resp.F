!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
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
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: do_resp_x_dir,&
                                             do_resp_y_dir,&
                                             do_resp_z_dir,&
                                             use_perd_none,&
                                             use_perd_xyz
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_wait
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE qs_collocate_density,            ONLY: calculate_rho_resp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_resp'

  PUBLIC :: resp_fit,&
            resp_type

  TYPE resp_type
   LOGICAL                                :: itc, rheavies, nonperiodic_sys
   INTEGER                                :: nres, ncons,&
                                             nrest_sec, ncons_sec,&
                                             npoints,stride(3),my_fit
   INTEGER, DIMENSION(:), POINTER         :: atom_surf_list
   INTEGER, DIMENSION(:,:), POINTER       :: fitpoints
   REAL(KIND=dp)                          :: rheavies_strength,&
                                             rmax,rmin,&
                                             length,eta
   REAL(KIND=dp),DIMENSION(3)             :: box_hi, box_low
   REAL(KIND=dp),DIMENSION(:),POINTER     :: range_surf 
   REAL(KIND=dp), DIMENSION(:),POINTER    :: rhs
   REAL(KIND=dp), DIMENSION(:, :),POINTER :: matrix
  END TYPE resp_type

CONTAINS

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
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(section_vals_type), POINTER         :: cons_section, input, &
                                                poisson_section, &
                                                resp_section, rest_section

    CALL timeset(routineN,handle)

    NULLIFY(logger,cell,subsys,particles,particle_set,input,resp_section,&
            cons_section,rest_section,poisson_section,resp_env)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, input=input, subsys=subsys,&
            particle_set=particle_set, cell=cell, error=error)
       resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                  error=error)
       CALL section_vals_get(resp_section, explicit=has_resp, error=error)
    END IF

    IF (.NOT. failure .AND. has_resp) THEN
       logger => cp_error_get_logger(error)
       poisson_section => section_vals_get_subs_vals(input,"DFT%POISSON",error=error)
       CALL section_vals_val_get(poisson_section,"PERIODIC",i_val=my_per,error=error)
       CALL create_resp_type(resp_env,error) 
       !initialize the RESP fitting, get all the keywords
       CALL init_resp(qs_env,resp_env,input,subsys,particle_set,cell,&
            resp_section,cons_section,rest_section,error)
       
       !print info
       CALL print_resp_parameter_info(qs_env,resp_env,my_per,error)
       
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
        CALL calc_resp_matrix_nonper(qs_env,resp_env,particles,cell,resp_env%matrix,&
                              resp_env%rhs,natom,error)   
       CASE(use_perd_xyz)
        CALL calc_resp_matrix_periodic(qs_env,resp_env,particles,cell,natom, error)
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
            natom,cons_section,particle_set,error)

       !solve system for the values of the charges and the lagrangian multipliers
       CALL DGETRF(nvar,nvar,resp_env%matrix,nvar,ipiv,info)
       CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

       CALL DGETRS('N',nvar,1,resp_env%matrix,nvar,ipiv,resp_env%rhs,nvar,info)
       CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

       CALL print_resp_charges(qs_env,resp_env,output_unit,natom,error)

       DEALLOCATE(ipiv, stat=stat)
       CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
       CALL resp_dealloc(resp_env,error)
       CALL cp_print_key_finished_output(output_unit,logger,resp_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE resp_fit
!****************************************************************************
!\brief creates the resp_type structure
!\param error variable to control error logging, stopping,...
!       see module cp_error_handling
!****************************************************************************
  SUBROUTINE create_resp_type(resp_env,error)
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_resp_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(resp_env)) CALL resp_dealloc(resp_env,error)
    ALLOCATE(resp_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

    IF (.NOT. failure) THEN
       NULLIFY(resp_env%atom_surf_list,&
               resp_env%range_surf,&
               resp_env%matrix,&
               resp_env%fitpoints,&
               resp_env%rhs)

    ENDIF

    resp_env%itc=.FALSE.
    resp_env%rheavies=.FALSE.
    resp_env%nonperiodic_sys=.FALSE.
 
    resp_env%box_hi=0.0_dp
    resp_env%box_low=0.0_dp

    resp_env%ncons=0
    resp_env%nres=0
    resp_env%npoints=0

  END SUBROUTINE create_resp_type
!****************************************************************************
!\brief deallocates the resp_type structure
!\param error variable to control error logging, stopping,...
!       see module cp_error_handling
!****************************************************************************
  SUBROUTINE resp_dealloc(resp_env, error)
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'resp_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(resp_env)) THEN
     IF (ASSOCIATED(resp_env%atom_surf_list)) THEN
      DEALLOCATE(resp_env%atom_surf_list, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     ENDIF
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
     DEALLOCATE(resp_env, stat=stat)
     CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    ENDIF
  END SUBROUTINE resp_dealloc
!****************************************************************************
!\brief intializes the resp fit. Getting the parameters
!****************************************************************************
  SUBROUTINE init_resp(qs_env,resp_env,input,subsys,particle_set,cell,resp_section,&
              cons_section,rest_section,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: resp_section, cons_section, &
                                                rest_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_resp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, i, j, &
                                                max_index, n_var, &
                                                num_surf_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: indexes, my_stride
    LOGICAL                                  :: explicit, failure, &
                                                index_in_range
    TYPE(section_vals_type), POINTER         :: nonperiodic_section, &
                                                periodic_section

    CALL timeset(routineN,handle)

    NULLIFY(indexes,my_stride,nonperiodic_section,periodic_section)
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
    CALL section_vals_get(nonperiodic_section, explicit=explicit, error=error)
    IF(explicit) THEN
      resp_env%nonperiodic_sys=.TRUE.
      CALL section_vals_val_get(nonperiodic_section,"RMIN",r_val=resp_env%rmin,&
                                                                    error=error)
      CALL section_vals_val_get(nonperiodic_section,"RMAX",r_val=resp_env%rmax,&
                                                                    error=error)
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

    ! get the parameter for periodic/surface systems
    CALL section_vals_get(periodic_section, explicit=explicit, error=error)
    IF(explicit) THEN
      CALL section_vals_val_get(periodic_section,"RANGE",r_vals=resp_env%range_surf,&
                                                                         error=error)
      CALL section_vals_val_get(periodic_section,"LENGTH", r_val=resp_env%length,&
                                                                      error=error)
      CALL section_vals_val_get(periodic_section,"SURF_DIRECTION",&
                                      i_val=resp_env%my_fit,error=error)
      IF(ANY(resp_env%range_surf<0.0_dp)) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Numbers in RANGE in PERIODIC_SYS cannot be negative.")
      ENDIF
      IF(resp_env%length<=EPSILON(0.0_dp)) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
             "Parameter LENGTH in PERIODIC_SYS has to be larger than zero.")
      ENDIF

      !list of atoms specifing the surface
      CALL section_vals_val_get(periodic_section,"ATOM_LIST",n_rep_val=n_var,error=error)
      num_surf_atom=0
      DO i=1,n_var
      CALL section_vals_val_get(periodic_section,"ATOM_LIST",i_rep_val=i,i_vals=indexes,&
                                                                             error=error)
      num_surf_atom=num_surf_atom + SIZE(indexes)
      ENDDO 
      ALLOCATE(resp_env%atom_surf_list(num_surf_atom),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      resp_env%atom_surf_list=0
      num_surf_atom=1
      DO i=1,n_var
       CALL section_vals_val_get(periodic_section,"ATOM_LIST",i_rep_val=i,i_vals=indexes,&
                                                                              error=error)
       resp_env%atom_surf_list(num_surf_atom:num_surf_atom+SIZE(indexes)-1)=indexes(:)
       num_surf_atom = num_surf_atom + SIZE(indexes)
      ENDDO
      !check atom list
      num_surf_atom=num_surf_atom-1
      max_index=SIZE(subsys%particles%els)
      CPPrecondition(SIZE(resp_env%atom_surf_list) /= 0,cp_failure_level,routineP,error,failure)
      index_in_range=(MAXVAL(resp_env%atom_surf_list)<= max_index)&
                           .AND.(MINVAL(resp_env%atom_surf_list) > 0)
      CPPostcondition(index_in_range,cp_failure_level,routineP,error,failure)
      DO i=1,num_surf_atom
       DO j=i+1,num_surf_atom
        atom_a=resp_env%atom_surf_list(i)
        atom_b=resp_env%atom_surf_list(j)
        IF(atom_a==atom_b) &
        CALL stop_program(routineN,moduleN,__LINE__,&
             "There are atoms doubled in surface atom list for RESP.")
       ENDDO
      ENDDO
    ENDIF

    ! get the parameters for the constraint and restraint sections    
    CALL section_vals_get(cons_section, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_get(cons_section,n_repetition=resp_env%ncons_sec,error=error)
    ELSE
       resp_env%ncons_sec=0
    ENDIF
    IF (explicit) THEN
       CALL section_vals_get(rest_section,n_repetition=resp_env%nrest_sec,error=error)
    ELSE
       resp_env%nrest_sec=0
    ENDIF
    resp_env%ncons=resp_env%ncons_sec
    resp_env%nres=resp_env%nrest_sec
    
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
 
!****************************************************************************
!\brief build matrix and vector for nonperiodic RESP fitting
!****************************************************************************
  SUBROUTINE calc_resp_matrix_nonper(qs_env,resp_env,particles,cell,matrix,&
                                                           rhs,natom,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rhs
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_resp_matrix_nonper', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bo(2,3), gbo(2,3), handle, i, &
                                                jx, jy, jz, k, l, m, np(3), &
                                                p, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dh(3,3), dvol, r(3), vec(3), &
                                                vec_pbc(3), vj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dist
    TYPE(pw_type), POINTER                   :: v_hartree_pw

    CALL timeset(routineN,handle)

    NULLIFY(v_hartree_pw)
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

    v_hartree_pw => qs_env%ks_env%v_hartree_rspace%pw
    bo=v_hartree_pw%pw_grid%bounds_local
    gbo=v_hartree_pw%pw_grid%bounds
    np=v_hartree_pw%pw_grid%npts
    dh=v_hartree_pw%pw_grid%dh
    dvol=v_hartree_pw%pw_grid%dvol
    ALLOCATE(dist(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

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
       DO i=1,natom
          vec=r-particles%els(i)%r
          vec_pbc(1) = vec(1) - cell%hmat(1,1)*ANINT(cell%h_inv(1,1)*vec(1))
          vec_pbc(2) = vec(2) - cell%hmat(2,2)*ANINT(cell%h_inv(2,2)*vec(2))
          vec_pbc(3) = vec(3) - cell%hmat(3,3)*ANINT(cell%h_inv(3,3)*vec(3))
          dist(i)=SQRT(SUM(vec_pbc**2))
       ENDDO
       ! check if the point is sufficiently close and sufficiently far. if OK, we can use
       ! the point for fitting, add/subtract 1.0E-13 to get rid of rounding errors when shifting atoms 
       IF (ALL(dist>resp_env%rmax-1.0E-13_dp).OR.ANY(dist<resp_env%rmin+1.0E-13_dp)) CYCLE
       resp_env%npoints = resp_env%npoints + 1

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

    CALL mp_sum(resp_env%npoints,v_hartree_pw%pw_grid%para%group)
    CALL mp_sum(matrix,v_hartree_pw%pw_grid%para%group)
    CALL mp_sum(rhs,v_hartree_pw%pw_grid%para%group)
    !weighted sum
    matrix=matrix/resp_env%npoints
    rhs=rhs/resp_env%npoints
    DEALLOCATE(dist,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)

    CALL timestop(handle)

  END SUBROUTINE calc_resp_matrix_nonper 
!****************************************************************************
!\brief build matrix and vector for periodic RESP fitting
!****************************************************************************
  SUBROUTINE calc_resp_matrix_periodic(qs_env,resp_env,particles,cell,natom,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_resp_matrix_periodic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ip, j, jx, jy, jz, &
                                                stat
    LOGICAL                                  :: failure
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
    CALL get_fitting_points(qs_env,resp_env,particles=particles,cell=cell,error=error)
    ALLOCATE(vpot(resp_env%npoints,natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)

    DO i=1,natom
     !collocate gaussian for each atom
      CALL pw_zero(rho_ga%pw, error=error)
      CALL calculate_rho_resp(rho_ga,qs_env,resp_env%eta,i,error)
     !calculate potential va and store the part needed for fitting in vpot
      CALL pw_zero(va_gspace%pw, error=error)
      CALL pw_poisson_solve(poisson_env,rho_ga%pw,vhartree=va_gspace%pw,error=error)
      CALL pw_zero(va_rspace%pw, error=error)
      CALL pw_transfer(va_gspace%pw,va_rspace%pw,error=error)
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

    CALL mp_sum(resp_env%npoints,para_env%group)
    CALL mp_sum(resp_env%matrix,para_env%group)
    CALL mp_sum(resp_env%rhs,para_env%group)
    !weighted sum
    resp_env%matrix=resp_env%matrix/resp_env%npoints
    resp_env%rhs=resp_env%rhs/resp_env%npoints

    CALL timestop(handle)

  END SUBROUTINE calc_resp_matrix_periodic

!****************************************************************************
!\brief get RESP fitting points for the periodic fitting
!****************************************************************************
  SUBROUTINE get_fitting_points(qs_env,resp_env,particles,cell,error)
 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_fitting_points', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2,3), gbo(2,3), handle, iatom, in_x, in_y, in_z, jx, jy, &
      jz, k, l, m, natom, now, np(3), output_unit, p, request, stat
    LOGICAL                                  :: failure, my_write
    LOGICAL, DIMENSION(:), POINTER           :: finishedloop
    REAL(KIND=dp)                            :: dh(3,3), dvol, r(3), &
                                                vec_pbc(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dist
    REAL, DIMENSION(:), POINTER              :: mepos_done
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_type), POINTER                   :: v_hartree_pw
    TYPE(section_vals_type), POINTER         :: input, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_write=.TRUE.
    NULLIFY(v_hartree_pw,para_env,logger,input,resp_section,finishedloop,mepos_done)   
 
    CALL get_qs_env(qs_env, input=input, para_env=para_env, error=error)

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
    np=v_hartree_pw%pw_grid%npts
    dh=v_hartree_pw%pw_grid%dh
    dvol=v_hartree_pw%pw_grid%dvol
    natom=SIZE(particles%els)

    IF(.NOT.ASSOCIATED(resp_env%fitpoints)) THEN
       now = 1000
       ALLOCATE(resp_env%fitpoints(3,now),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    ELSE
       now = SIZE(resp_env%fitpoints,2)
    END IF

    ALLOCATE(finishedloop(para_env%num_pe),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    finishedloop=.FALSE.
    ALLOCATE(dist(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    dist=0.0_dp

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
         DO m=1,natom
            vec_pbc = pbc(r,particles%els(m)%r,cell)
            dist(m)=SQRT(SUM(vec_pbc**2))
         ENDDO
         IF (ALL(dist>resp_env%rmax-1.0E-13_dp).OR.ANY(dist<resp_env%rmin+1.0E-13_dp)) CYCLE
       ELSE
         DO m=1,SIZE(resp_env%atom_surf_list)
            in_z=0
            in_y=0
            in_x=0
            iatom=resp_env%atom_surf_list(m)
            vec_pbc = pbc(particles%els(iatom)%r,r,cell)
            SELECT CASE(resp_env%my_fit)
            !subtract 1.0E-13 to get rid of rounding errors when shifting atoms 
            CASE(do_resp_x_dir)
             IF(ABS(vec_pbc(3))<resp_env%length-1.0E-13_dp)       in_z=1
             IF(ABS(vec_pbc(2))<resp_env%length-1.0E-13_dp)       in_y=1
             IF(vec_pbc(1)>resp_env%range_surf(1)+1.0E-13_dp.AND.&              
                 vec_pbc(1)<resp_env%range_surf(2)-1.0E-13_dp)    in_x=1 
            CASE(do_resp_y_dir)
             IF(ABS(vec_pbc(3))<resp_env%length-1.0E-13_dp)       in_z=1
             IF(vec_pbc(2)>resp_env%range_surf(1)+1.0E-13_dp.AND.&              
                 vec_pbc(2)<resp_env%range_surf(2)-1.0E-13_dp)    in_y=1 
             IF(ABS(vec_pbc(1))<resp_env%length-1.0E-13_dp)       in_x=1
            CASE(do_resp_z_dir)
             IF(vec_pbc(3)>resp_env%range_surf(1)+1.0E-13_dp.AND.&              
                 vec_pbc(3)<resp_env%range_surf(2)-1.0E-13_dp)    in_z=1 
             IF(ABS(vec_pbc(2))<resp_env%length-1.0E-13_dp)       in_y=1
             IF(ABS(vec_pbc(1))<resp_env%length-1.0E-13_dp)       in_x=1
            END SELECT
            IF(in_z*in_y*in_x==1) EXIT
         ENDDO
         IF (BTEST(cp_print_key_should_output(logger%iter_info,&
                   resp_section,"PRINT%COORD_FIT_POINTS",error=error),&
                   cp_p_file))THEN
           CALL print_fitting_points(qs_env,resp_env,r,in_x,in_y,in_z,&
                output_unit,my_write,finishedloop,error)
         ENDIF
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
 
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
            resp_section,"PRINT%COORD_FIT_POINTS",error=error),&
            cp_p_file))THEN
     !procs must tell proc 0 that they finished loops
     IF(para_env%mepos/=0) THEN
      ALLOCATE(mepos_done(3))
      CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
      mepos_done=-2.0_dp ! -2.0 is arbitrary, acts as a flag
      CALL mp_isend(msgin=mepos_done,dest=0,comm=para_env%group,request=request)
      CALL mp_wait(request)
      DEALLOCATE(mepos_done)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
     !if proc 0 finishes loops earlier than other procs, it has to recv their data here 
     ELSE
      finishedloop(1)=.TRUE. !proc 0 finished
      DO
       IF(.NOT.ALL(finishedloop.eqv..TRUE.)) THEN
        CALL print_fitting_points(qs_env,resp_env,r,0,0,0,&
             output_unit,my_write,finishedloop,error)
       ELSE
        EXIT
       ENDIF
      ENDDO
     ENDIF
    ENDIF
    CALL cp_print_key_finished_output(output_unit,logger,resp_section,&
                       "PRINT%COORD_FIT_POINTS", error=error)
    DEALLOCATE(finishedloop,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    DEALLOCATE(dist,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)

    CALL timestop(handle)
 
  END SUBROUTINE get_fitting_points

!****************************************************************************
!\brief calculate vector rhs
!****************************************************************************
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
!****************************************************************************
!\brief print the atom coordinates and the coordinates of the fitting points
!       to an xyz file
!****************************************************************************
  SUBROUTINE print_fitting_points(qs_env,resp_env,r,in_x,in_y,in_z,output_unit,&
                                  my_write,finishedloop,error)
  
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    INTEGER, INTENT(IN)                      :: in_x, in_y, in_z, output_unit
    LOGICAL, INTENT(INOUT)                   :: my_write
    LOGICAL, DIMENSION(:), POINTER           :: finishedloop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_fitting_points', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, i, iatom, m, req(2), &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: conv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_proc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tmp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: input, print_key, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(para_env,input,logger,resp_section,print_key,particle_set,tmp)   
 
    CALL get_qs_env(qs_env, input=input, para_env=para_env,&
                       particle_set=particle_set,error=error)
 
    resp_section => section_vals_get_subs_vals(input,"PROPERTIES%RESP",&
                                                error=error)
    print_key => section_vals_get_subs_vals(resp_section,&
                                            "PRINT%COORD_FIT_POINTS",&
                                                error=error)
    logger => cp_error_get_logger(error)
    conv=cp_unit_from_cp2k(1.0_dp,"angstrom",error=error)
 
    ALLOCATE(r_proc(para_env%num_pe,3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    !every r is >= 0 by definition, print only if r_proc >= 0 in the end
    r_proc=-1.0_dp

    !this condition is NOT true for every proc 
    IF(in_z*in_y*in_x==1) THEN
      r_proc(para_env%mepos+1,:)=r(:)*conv
    ENDIF

    !sending it all to proc 0 making the printing
    IF(para_env%mepos==0) THEN
      ALLOCATE(tmp(3),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
      DO i=2,para_env%num_pe
       IF(finishedloop(i)) CYCLE
       tmp(:)=r_proc(i,:)
       CALL mp_irecv(msgout=tmp,source=i-1,comm=para_env%group,request=req(1))
       CALL mp_wait(req(1))
       !tmp=-2.0 (mepos_done) is received when proc i finished the 3 loops
       !tmp=-1.0 is received if in_z*in_y*in_x=0, otherwise tmp >= 0
       IF(ALL(tmp >= -1.0_dp)) THEN
          r_proc(i,:)=tmp
       ELSE
          finishedloop(i) = .TRUE.
       ENDIF
      ENDDO
    ELSE
     ALLOCATE(tmp(3),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
     tmp(:)=r_proc(para_env%mepos+1,:)
     CALL mp_isend(msgin=tmp,dest=0,comm=para_env%group,request=req(2))
     CALL mp_wait(req(2))
    ENDIF

    IF(ASSOCIATED(tmp)) THEN
      DEALLOCATE(tmp, stat=stat)
      CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    ENDIF    
       
    IF(output_unit>0) THEN
     filename=cp_print_key_generate_filename(logger,&
              print_key, extension=".xyz",&
              my_local=.FALSE.,error=error)
     IF(my_write) THEN
      WRITE(unit=output_unit,FMT="(I12,A,/)") SIZE(particle_set),' + nr fit points'
      DO iatom=1,SIZE(particle_set)
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                            element_symbol=element_symbol)
       WRITE(UNIT=output_unit,FMT="(A,1X,3F10.5)") element_symbol,&
                                               particle_set(iatom)%r(1:3)*conv
       my_write=.FALSE.  
      ENDDO
     ENDIF     
     DO m=1,para_env%num_pe
      IF(ANY(r_proc(m,:)>=0.0_dp)) THEN 
       WRITE(UNIT=output_unit,FMT="(A,2X,3F10.5)") "X", r_proc(m,1), r_proc(m,2), r_proc(m,3)
      ENDIF
     ENDDO
    ENDIF
 
    DEALLOCATE(r_proc,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)

    CALL timestop(handle)
 
  END SUBROUTINE print_fitting_points
!****************************************************************************
!\brief print input information
!****************************************************************************
  SUBROUTINE print_resp_parameter_info(qs_env,resp_env,my_per,error)
   
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    INTEGER, INTENT(IN)                      :: my_per
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_resp_parameter_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
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
     WRITE(output_unit,'(T3,A,T75,I6)') "Number of explicit constraints: ",resp_env%ncons_sec
     WRITE(output_unit,'(T3,A,T75,I6)') "Number of explicit restraints: ",resp_env%nrest_sec
     WRITE(output_unit,'(T3,A,T80,A)')  "Constrain total charge ",MERGE("T","F",resp_env%itc)
     WRITE(output_unit,'(T3,A,T80,A)')  "Restrain heavy atoms ",MERGE("T","F",resp_env%rheavies)
     IF (resp_env%rheavies) THEN
        WRITE(output_unit,'(T3,A,T71,F10.5)') "Heavy atom restraint strength: ",&
                                                                    resp_env%rheavies_strength
     ENDIF
     WRITE(output_unit,'(T3,A,T66,3I5)') "Stride: ",resp_env%stride
     IF(resp_env%nonperiodic_sys) THEN
       WRITE(output_unit,'(T3,A,T71,F10.5)')  "Rmin [angstrom]: " ,resp_env%rmin*conv
       WRITE(output_unit,'(T3,A,T71F10.5)')  "Rmax [angstrom]: ",resp_env%rmax*conv
       WRITE(output_unit,'(T3,A,T51,3F10.5)') "Box min [angstrom]: ",resp_env%box_low(1:3)*conv
       WRITE(output_unit,'(T3,A,T51,3F10.5)') "Box max [angstrom]: ",resp_env%box_hi(1:3)*conv
     ELSE
       WRITE(output_unit,'(2X,A,F10.5)')  "Index of atoms defining the surface: "
       WRITE(output_unit,'(7X,10I6)')  resp_env%atom_surf_list
       WRITE(output_unit,'(T3,A,T61,2F10.5)')&
                          "Range for sampling above the surface [angstrom]:",&
                           resp_env%range_surf(1:2)*conv
       WRITE(output_unit,'(T3,A,T71,F10.5)')  "Length of sampling box above each"//&
                                              " surface atom [angstrom]: ",resp_env%length*conv
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
!****************************************************************************
!\brief add restraints and constraints
!****************************************************************************
  SUBROUTINE add_restraints_and_constraints(qs_env,resp_env,rest_section,&
                                        natom,cons_section,particle_set,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(resp_type), POINTER                 :: resp_env
    TYPE(section_vals_type), POINTER         :: rest_section
    INTEGER, INTENT(IN)                      :: natom
    TYPE(section_vals_type), POINTER         :: cons_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'add_restraints_and_constraints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, i, j, &
                                                k, l, m, ncons_i, ncons_v, z
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, index_in_range
    REAL(KIND=dp)                            :: strength, TARGET
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atom_coef

    CALL timeset(routineN,handle)

    NULLIFY(atom_coef,atom_list)
    failure=.FALSE.
 
    ! add the restraints
    DO i=1,resp_env%nrest_sec
       CALL section_vals_val_get(rest_section,"ATOM_LIST",i_rep_section=i,i_vals=atom_list,error=error)
       CALL section_vals_val_get(rest_section,"ATOM_COEF",i_rep_section=i,r_vals=atom_coef,error=error)
       CALL section_vals_val_get(rest_section,"TARGET",i_rep_section=i,r_val=TARGET,error=error)
       CALL section_vals_val_get(rest_section,"STRENGTH",i_rep_section=i,r_val=strength,error=error)
       CPPrecondition(SIZE(atom_list)==SIZE(atom_coef),cp_failure_level,routineP,error,failure)
       !check atom_list
       index_in_range=(MAXVAL(atom_list) <= natom).AND.(MINVAL(atom_list) > 0)
       CPPostcondition(index_in_range,cp_failure_level,routineP,error,failure)
       DO l=1,SIZE(atom_list)
        DO j=l+1,SIZE(atom_list)
         atom_a=atom_list(l)
         atom_b=atom_list(j)
         IF(atom_a==atom_b) &
         CALL stop_program(routineN,moduleN,__LINE__,&
              "There are atoms doubled in atom list of restraint section.")
        ENDDO
       ENDDO
       DO m=1,SIZE(atom_list)
          DO k=1,SIZE(atom_list)
             resp_env%matrix(atom_list(m),atom_list(k))=resp_env%matrix(atom_list(m),atom_list(k))+ &
                                                atom_coef(m)*atom_coef(k)*2.0_dp*strength
          ENDDO
          resp_env%rhs(atom_list(m))=resp_env%rhs(atom_list(m))+2.0_dp*TARGET*strength*atom_coef(m)
       ENDDO
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
    ncons_i=0
    IF (resp_env%itc) THEN
       ncons_i=ncons_i+1
       ncons_v=natom+ncons_i
       resp_env%matrix(1:natom,ncons_v)=1.0_dp
       resp_env%matrix(ncons_v,1:natom)=1.0_dp
       resp_env%rhs(ncons_v)=qs_env%dft_control%charge
    ENDIF

    DO i=1,resp_env%ncons_sec
       ncons_i=ncons_i+1
       ncons_v=ncons_v+1
       CALL section_vals_val_get(cons_section,"ATOM_LIST",i_rep_section=i,i_vals=atom_list,error=error)
       CALL section_vals_val_get(cons_section,"ATOM_COEF",i_rep_section=i,r_vals=atom_coef,error=error)
       CALL section_vals_val_get(cons_section,"TARGET",i_rep_section=i,r_val=TARGET,error=error)
       CPPrecondition(SIZE(atom_list)==SIZE(atom_coef),cp_failure_level,routineP,error,failure)
       !check atom_list
       index_in_range=(MAXVAL(atom_list) <= natom).AND.(MINVAL(atom_list) > 0)
       CPPostcondition(index_in_range,cp_failure_level,routineP,error,failure)
       DO l=1,SIZE(atom_list)
        DO j=l+1,SIZE(atom_list)
         atom_a=atom_list(l)
         atom_b=atom_list(j)
         IF(atom_a==atom_b) &
         CALL stop_program(routineN,moduleN,__LINE__,&
              "There are atoms doubled in atom list of constraint section.")
        ENDDO
       ENDDO
       DO m=1,SIZE(atom_list)
          resp_env%matrix(atom_list(m),ncons_v)=atom_coef(m)
          resp_env%matrix(ncons_v,atom_list(m))=atom_coef(m)
       ENDDO
       resp_env%rhs(ncons_v)=TARGET
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE add_restraints_and_constraints
!****************************************************************************
!\brief print RESP charges to an extra file or to the normal output file
!****************************************************************************
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
    TYPE(section_vals_type), POINTER         :: input, print_key, resp_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(particle_set,input,logger,resp_section,print_key)
 
    CALL get_qs_env(qs_env,input=input,particle_set=particle_set,error=error)

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
     CALL print_atomic_charges(particle_set,output_file,title="RESP charges:",&
                                             atomic_charges=resp_env%rhs(1:natom))
     IF(output_runinfo>0) WRITE(output_runinfo,'(2X,A,/)')  "PRINTED RESP CHARGES TO FILE"
     ENDIF
 
     CALL cp_print_key_finished_output(output_file,logger,resp_section,&
                       "PRINT%RESP_CHARGES_TO_FILE", error=error)
    ELSE
     CALL print_atomic_charges(particle_set,output_runinfo,title="RESP charges:",&
                                             atomic_charges=resp_env%rhs(1:natom))
    ENDIF

    CALL timestop(handle)
  END SUBROUTINE print_resp_charges
END MODULE qs_resp
