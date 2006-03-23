!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/pint_types *
!!
!!   NAME
!!     pint_types
!!
!!   FUNCTION
!!     types to performs a path integral run
!!
!!   NOTES
!!     quick & dirty rewrite of my python program
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2006 created
!!
!!   SOURCE
!****************************************************************************
MODULE pint_types
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE input_constants,                 ONLY: dump_dcd
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             next_rng_seed,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE replica_methods,                 ONLY: rep_env_calc_e_f,&
                                             rep_env_create
  USE replica_types,                   ONLY: rep_env_release,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_types'
  
  INTEGER, SAVE, PRIVATE :: last_pint_id=0, last_staging_id=0

  PUBLIC :: pint_env_type, staging_env_type
  PUBLIC :: pint_env_create, pint_retain, pint_release,&
       pint_x2u, pint_u2x, pint_f2uf, pint_calc_uf_h,&
       pint_calc_f, pint_init_v, pint_write_centroids, pint_write_ener,&
       pint_write_restart, pint_step, pint_start, pint_run
  PUBLIC :: staging_env_create, staging_x2u, staging_u2x, staging_f2uf,&
       staging_init_masses
!***
!****************************************************************************

  !!****s* pint_types/pint_env_type *
  !!
  !!   NAME
  !!     pint_env_type
  !!
  !!   FUNCTION
  !!     environment for a path integral run
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE pint_env_type
     INTEGER :: ref_count, id_nr, p,nnos,nrespa, iter,ndim
     REAL(kind=dp) :: v_tol,kT,beta,dt,&
          e_pot_h, e_kin_beads, e_pot_t, e_kin_t,t
     TYPE(replica_env_type), POINTER :: replicas
     TYPE(section_vals_type), POINTER :: input
     TYPE(staging_env_type), POINTER :: staging_env
     TYPE(rng_stream_type), POINTER :: randomG
     REAL(kind=dp), DIMENSION(:), POINTER :: mass,e_pot_bead
     REAL(kind=dp), DIMENSION(:,:), POINTER :: x,v,f,mass_beads,&
          mass_fict,ux,uv,uv_t,uv_new,uf,uf_h
     REAL(kind=dp), DIMENSION(:,:,:), POINTER :: tx,tv,tv_t,tv_old,tv_new,tf
     REAL(kind=dp), DIMENSION(:), POINTER :: Q ! dim p, make it (p,ndim)?
  END TYPE pint_env_type
  !!***
  !****************************************************************************

  !!****s* pint_types/staging_env_type *
  !!
  !!   NAME
  !!     staging_env_type
  !!
  !!   FUNCTION
  !!     data to perform the staging transformation
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE staging_env_type
     INTEGER :: id_nr,ref_count,nseg,j,p
     REAL(kind=dp) w_p,w_j,Q_stage,Q_end
  END TYPE staging_env_type
  !!***
  !****************************************************************************

CONTAINS

!!****f* pint_types/pint_env_create *
!!
!!   NAME
!!     pint_env_create
!!
!!   FUNCTION
!!     creates a path integral envirnment
!!
!!   NOTES
!!     might return an unassociated pointer in parallel on the processors
!!     that are not needed
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE pint_env_create(pint_env,input,para_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, ib, idim, idir, ierr, &
                                                nrep, prep, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: mass
    REAL(kind=dp), DIMENSION(3, 2)           :: seed
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER         :: pint_section, staging_section

    failure=.FALSE.
    NULLIFY(f_env,subsys,particles)

    CPPrecondition(.NOT.ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
    CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(rep_env)
       CALL section_vals_val_get(pint_section,"p",i_val=nrep,error=error)
       CALL rep_env_create(rep_env, para_env=para_env, input=input,&
            nrep=nrep,prep=prep,error=error)
       IF (.NOT. ASSOCIATED(rep_env)) RETURN

       ALLOCATE(pint_env,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       last_pint_id=last_pint_id+1
       pint_env%id_nr=last_pint_id
       pint_env%ref_count=1
       NULLIFY(pint_env%replicas,pint_env%input,pint_env%staging_env)
       pint_section => section_vals_get_subs_vals(input,"MOTION%PINT",&
            error=error)
       CALL section_vals_val_get(pint_section,"p",i_val=pint_env%p,error=error)
       CALL section_vals_val_get(pint_section,"proc_per_replica",&
            i_val=prep,error=error)
       pint_env%replicas => rep_env
       pint_env%ndim=rep_env%ndim
       CALL section_vals_val_get(staging_section,"nnos",i_val=pint_env%nnos,&
            error=error)
       CALL section_vals_val_get(staging_section,"nrespa",i_val=pint_env%nrespa,&
            error=error)
       CALL section_vals_val_get(staging_section,"Temp",r_val=pint_env%kT,&
            error=error)
       CALL section_vals_val_get(staging_section,"dt",r_val=pint_env%dt,&
            error=error)
       CALL section_vals_val_get(staging_section,"v_tol_nose",r_val=pint_env%v_tol,&
            error=error)
       pint_env%beta=1._dp/pint_env%kT

       seed(:,:) = next_rng_seed()
       CALL create_rng_stream(pint_env%randomG,&
                              name="pint_randomG",&
                              distribution_type=GAUSSIAN,&
                              extended_precision=.TRUE.,&
                              seed=seed)

       ALLOCATE(pint_env%e_pot_bead(pint_env%p),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       pint_env%e_pot_bead=0._dp
       pint_env%e_pot_h=0._dp
       pint_env%e_kin_beads=0._dp
       pint_env%e_pot_t=0._dp
       pint_env%e_kin_t=0._dp
       pint_env%iter=0
       pint_env%t=0._dp

       ALLOCATE(&
            pint_env%x(pint_env%p,pint_env%ndim),          &
            pint_env%v(pint_env%p,pint_env%ndim),          &
            pint_env%f(pint_env%p,pint_env%ndim),          &
            pint_env%ux(pint_env%p,pint_env%ndim),         &
            pint_env%uv(pint_env%p,pint_env%ndim),         &
            pint_env%uv_t(pint_env%p,pint_env%ndim),       &
            pint_env%uv_new(pint_env%p,pint_env%ndim),     &
            pint_env%uf(pint_env%p,pint_env%ndim),         &
            pint_env%uf_h(pint_env%p,pint_env%ndim),       &
            pint_env%tx(pint_env%nnos,pint_env%p,pint_env%ndim),    &
            pint_env%tv(pint_env%nnos,pint_env%p,pint_env%ndim),    &
            pint_env%tv_t(pint_env%nnos,pint_env%p,pint_env%ndim),  &
            pint_env%tv_old(pint_env%nnos,pint_env%p,pint_env%ndim),&
            pint_env%tv_new(pint_env%nnos,pint_env%p,pint_env%ndim),& 
            pint_env%tf(pint_env%nnos,pint_env%p,pint_env%ndim), stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       DO ib=1,pint_env%p
          DO idim=1,pint_env%ndim
             pint_env%x(ib,idim)=pint_env%replicas%r(idim,ib)
          END DO
       END DO
       pint_env%v        = 0._dp
       pint_env%f        = 0._dp
       pint_env%ux       = 0._dp
       pint_env%uv       = 0._dp
       pint_env%uv_t     = 0._dp
       pint_env%uv_new   = 0._dp
       pint_env%uf       = 0._dp
       pint_env%uf_h     = 0._dp
       pint_env%tx       = 0._dp
       pint_env%tv       = 0._dp
       pint_env%tv_t     = 0._dp
       pint_env%tv_old   = 0._dp
       pint_env%tv_new   = 0._dp
       pint_env%tf       = 0._dp

       staging_section => section_vals_get_subs_vals(input,"MOTION%PINT%STAGING",&
            error=error)
       CALL staging_env_create(pint_env%staging_env,staging_section,p=pint_env%p,&
            ndim=pint_env%ndim,kT=pint_env%kT, error=error)
       ALLOCATE(pint_env%mass(pint_env%ndim),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
            f_env=f_env,new_error=new_error, failure=failure)
       logger => cp_error_get_logger(new_error)
       CALL force_env_get(force_env=f_env%force_env,&
            subsys=subsys,&
            error=new_error)
       CALL cp_subsys_get(subsys(1)%subsys,&
            particles=particles,&
            error=new_error)

       idim=0
       DO iat=1,pint_env%ndim
          CALL get_atomic_kind(particles%els(iat)%atomic_kind,mass=mass)
          DO idir=1,3
             idim=idim+1
             pint_env%mass(idim)=mass
          END DO
       END DO
       CALL cp_assert(logger%iter_info%iteration(2)==pint_env%iter,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "md & force_eval lost sychro "//&
            CPSourceFileRef,&
            error,failure)
       CALL f_env_rm_defaults(f_env,new_error,ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

       ALLOCATE(pint_env%Q(pint_env%p),&
            pint_env%mass_beads(pint_env%p,pint_env%ndim),&
            pint_env%mass_fict(pint_env%p,pint_env%ndim),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL staging_init_masses(pint_env%staging_env,mass=pint_env%mass,&
            mass_beads=pint_env%mass_beads,mass_fict=pint_env%mass_fict,&
            Q=pint_env%Q,error=error)

    END IF
  END SUBROUTINE pint_env_create
!***************************************************************************

!!****f* pint_types/pint_retain [1.0] *
!!
!!   NAME
!!     pint_retain
!!
!!   FUNCTION
!!     retains a path integral environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the pint_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE pint_retain(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
     pint_env%ref_count=pint_env%ref_count+1
  END IF
END SUBROUTINE pint_retain
!***************************************************************************

!!****f* pint_types/pint_release [1.0] *
!!
!!   NAME
!!     pint_release
!!
!!   FUNCTION
!!     releases a path integral environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the pint_env to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE pint_release(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(pint_env)) THEN
     CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
     pint_env%ref_count=pint_env%ref_count-1
     IF (pint_env%ref_count==0) THEN
        CALL rep_env_release(pint_env%replicas,error=error)
        CALL section_vals_release(pint_env%input,error=error)
        CALL staging_release(pint_env%staging_env,error=error)
        CALL delete_rng_stream(pint_env%randomG,error=error)

        DEALLOCATE(pint_env%mass,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%e_pot_bead,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DEALLOCATE(pint_env%x,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%v,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%f,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%mass_beads,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%mass_fict,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%ux,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uv,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uv_t,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uf,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uf_h,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DEALLOCATE(pint_env%tx,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%tv,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%tv_t,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%tv_old,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%tv_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%tf,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DEALLOCATE(pint_env%Q,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DEALLOCATE(pint_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(pint_env)
END SUBROUTINE pint_release
!***************************************************************************

!!****f* pint_types/staging_retain [1.0] *
!!
!!   NAME
!!     staging_retain
!!
!!   FUNCTION
!!     retains a staging_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the staging_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE staging_retain(staging_env,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(staging_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(staging_env%ref_count>0,cp_failure_level,routineP,error,failure)
     staging_env%ref_count=staging_env%ref_count+1
  END IF
END SUBROUTINE staging_retain
!***************************************************************************

!!****f* pint_types/staging_release [1.0] *
!!
!!   NAME
!!     staging_release
!!
!!   FUNCTION
!!     releases the staging environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the staging_env to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE staging_release(staging_env,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(staging_env)) THEN
     CPPrecondition(staging_env%ref_count>0,cp_failure_level,routineP,error,failure)
     staging_env%ref_count=staging_env%ref_count-1
     IF (staging_env%ref_count==0) THEN
        DEALLOCATE(staging_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(staging_env)
END SUBROUTINE staging_release
!***************************************************************************

!!****f* pint_types/staging_init_masses [1.0] *
!!
!!   NAME
!!     staging_init_masses
!!
!!   FUNCTION
!!     initializes the masses and fictitious masses compatibly with the
!!     staging information
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the definition of the staging transformation
!!     - masses: *input* the masses of the particles
!!     - mass_beads: masses of the beads
!!     - mass_fict: the fictitious masses
!!     - Q: masses of the nose thermostats
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE staging_init_masses(staging_env,mass,mass_beads,mass_fict,&
     Q,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    REAL(kind=dp), DIMENSION(:), INTENT(in)  :: mass
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL                  :: mass_beads, mass_fict
    REAL(kind=dp), DIMENSION(:), &
      INTENT(out), OPTIONAL                  :: Q
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_init_masses', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, ib, iseg, stat
    LOGICAL                                  :: failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: scal

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     IF (PRESENT(Q)) THEN
        Q=staging_env%Q_stage
        DO i=1,staging_env%p,staging_env%j
           Q(i)=staging_env%Q_end
        END DO
     END IF
     IF (PRESENT(mass_beads).OR.PRESENT(mass_fict)) THEN
        
        ALLOCATE(scal(staging_env%p),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        DO iseg=1,staging_env%nseg
           DO i=1,staging_env%j ! check order!!!
              scal(staging_env%j*(iseg-1)+i)=REAL(i,dp)/REAL(MAX(1,i-1),dp)
           END DO
        END DO
     !   scal=zeros(staging_env%j,Float64)
     !   divide(arange(2,staging_env%j+1,typecode=Float64),
     !          arange(1,staging_env%j,typecode=Float64),scal[1:])
     !   scal[0]=1.
     !   scal=outerproduct(ones(staging_env%nseg),scal)

        IF (PRESENT(mass_beads)) THEN
           DO iat=1,SIZE(mass)
              DO ib=1,staging_env%p
                 mass_beads(ib,iat)=scal(ib)*mass(iat)
              END DO
           END DO
        END IF
        IF (PRESENT(mass_fict)) THEN
           DO iat=1,SIZE(mass)
              DO ib=1,staging_env%p
                 mass_fict(ib,iat)=scal(ib)*mass(iat)
              END DO
           END DO
        END IF
     END IF
  END IF
END SUBROUTINE staging_init_masses
!***************************************************************************

!!****f* pint_types/staging_env_create *
!!
!!   NAME
!!     staging_env_create
!!
!!   FUNCTION
!!     creates the data needed for a staging transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE staging_env_create(staging_env,staging_section,p,ndim,kT,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    TYPE(section_vals_type), POINTER         :: staging_section
    INTEGER, INTENT(in)                      :: p, ndim
    REAL(kind=dp), INTENT(in)                :: kT
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  NULLIFY(staging_section)

  CPPrecondition(.NOT.ASSOCIATED(staging_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     ALLOCATE(staging_env,stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     last_staging_id=last_staging_id+1
     staging_env%id_nr=last_staging_id
     staging_env%ref_count=1

     CALL section_vals_val_get(staging_section,"j",i_val=staging_env%j,&
          error=error)
     CALL section_vals_val_get(staging_section,"Q_end",i_val=staging_env%j,&
          error=error)
     staging_env%p=p
     staging_env%nseg=staging_env%p/staging_env%j
     
     staging_env%w_p=SQRT(REAL(staging_env%p,dp))*kT
     staging_env%w_j=kT*SQRT(REAL(staging_env%nseg,dp))
     staging_env%Q_stage=kT/staging_env%w_p**2
     IF (staging_env%Q_end<=0._dp) THEN
        staging_env%Q_end=staging_env%j*staging_env%Q_stage
     END IF
  END IF
END SUBROUTINE staging_env_create
!***************************************************************************

!!****f* pint_types/pint_x2u *
!!
!!   NAME
!!     pint_x2u
!!
!!   FUNCTION
!!     Transforms from the x into the u variables
!!     (at the moment a staging transformation for the positions)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the path integral environment
!!     - ux: will contain the u variable (defaults to pint_env%ux)
!!     - x: the positions to transform (defaults to pint_env%x)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_x2u(pint_env,ux,x,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: x
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_x2u', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_ux, my_x

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_x => pint_env%x
     my_ux => pint_env%ux
     IF (PRESENT(x)) my_x => x
     IF (PRESENT(ux)) my_ux => ux
     CPPrecondition(ASSOCIATED(my_ux),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(my_x),cp_failure_level,routineP,error,failure)
     
     CALL staging_x2u(pint_env%staging_env,ux=my_ux,x=my_x,error=error)
  END IF
END SUBROUTINE pint_x2u
!***************************************************************************

!!****f* pint_types/staging_x2u *
!!
!!   NAME
!!     staging_x2u
!!
!!   FUNCTION
!!     Transforms from the x into the u variables using a staging
!!     transformation for the positions
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the environment for the staging transformation
!!     - ux: will contain the u variable
!!     - x: the positions to transform
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE staging_x2u(staging_env,ux,x,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: x
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_x2u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, s
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(staging_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(staging_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ux=x
     DO s=0,staging_env%nseg-1
        DO k=2,staging_env%j
           ux(staging_env%j*s+k,:)=ux(staging_env%j*s+k,:)&
                -((REAL(k-1,dp)/REAL(k,dp)&
                *x(MODULO((staging_env%j*s+k+1),staging_env%p),:)+ &
                x(staging_env%j*s+1,:)/REAL(k,dp)))
        END DO
     END DO
  END IF
END SUBROUTINE staging_x2u
!***************************************************************************

!!****f* pint_types/pint_u2x *
!!
!!   NAME
!!     pint_u2x
!!
!!   FUNCTION
!!     transform from the u variable to the x (inverse of x2u)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: path integral environment
!!     - ux: the u variable (positions to be backtransformed)
!!     - x: will contain the positions
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_u2x(pint_env,ux,x,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: x
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_u2x', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_ux, my_x

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_x => pint_env%x
     my_ux => pint_env%ux
     IF (PRESENT(x)) my_x => x
     IF (PRESENT(ux)) my_ux => ux
     CPPrecondition(ASSOCIATED(my_ux),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(my_x),cp_failure_level,routineP,error,failure)
     
     CALL staging_u2x(pint_env%staging_env,ux=my_ux,x=my_x,error=error)
  END IF
END SUBROUTINE pint_u2x
!***************************************************************************

!!****f* pint_types/staging_u2x *
!!
!!   NAME
!!     staging_u2x
!!
!!   FUNCTION
!!     transform from the u variable to the x (back staging transformation
!!     for the positions)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the environment for the staging transformation
!!     - ux: the u variable (positions to be backtransformed)
!!     - x: will contain the positions
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE staging_u2x(staging_env,ux,x,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: x
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_u2x', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ist, j, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii, jjj
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: const, const2

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(staging_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(staging_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     j=staging_env%j
     const=REAL(j-1,dp)/REAL(j,dp)
     const2=1._dp/REAL(j,dp)
     ALLOCATE(iii(staging_env%nseg),jjj(staging_env%nseg),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     DO i=1,staging_env%nseg
        iii(i)=staging_env%j*(i-1)+1 !first el
     END DO
     DO i=1,staging_env%nseg-1
        jjj(i)=iii(i)+j ! next first el (pbc)
     END DO
     jjj(staging_env%nseg)=1
            
     x=ux
     DO i=1,staging_env%nseg
        x(j-1+iii(i),:)=x(j-1+iii(i),:)+&
             const*ux(jjj(i),:)+ux(iii(i),:)*const2
     END DO
     DO ist=1,staging_env%nseg
        DO i=staging_env%j-2,2,-1
           x(i+iii(ist),:)=x(i+iii(ist),:)+&
                REAL(i-1,dp)/REAL(i,dp)*x(i+iii(ist)+1,:)&
                +ux(iii(ist),:)/REAL(i,dp)
        END DO
     END DO
  END IF
END SUBROUTINE staging_u2x
!***************************************************************************

!!****f* pint_types/pint_f2uf *
!!
!!   NAME
!!     pint_f2uf
!!
!!   FUNCTION
!!     transformation x to u for the forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the path integral environment
!!     - uf: will contain the forces after for the tranformed variable
!!     - f: the forces to transform
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_f2uf(pint_env,uf,f,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: uf
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: f
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_f2uf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_f, my_uf

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_f => pint_env%f
     my_uf => pint_env%uf
     IF (PRESENT(f)) my_f => f
     IF (PRESENT(uf)) my_uf => uf
     CPPrecondition(ASSOCIATED(my_uf),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(my_f),cp_failure_level,routineP,error,failure)
     
     CALL staging_f2uf(pint_env%staging_env,uf=my_uf,f=my_f,error=error)
     my_uf=my_uf/pint_env%mass_fict
  END IF
END SUBROUTINE pint_f2uf
!***************************************************************************

!!****f* pint_types/staging_f2uf *
!!
!!   NAME
!!     staging_f2uf
!!
!!   FUNCTION
!!     staging transformation for the forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the environment for the staging transformation
!!     - uf: will contain the forces after for the tranformed variable
!!     - f: the forces to transform
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE staging_f2uf(staging_env,uf,f,error)
    TYPE(staging_env_type), POINTER          :: staging_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: uf
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: f
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_f2uf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idim, ij, ist, k, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii, jjj, kkk
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: const, sum_f

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(staging_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(staging_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     const=REAL(staging_env%j-1,dp)/REAL(staging_env%j,dp)
     ALLOCATE(iii(staging_env%j),jjj(staging_env%j),&
          kkk(staging_env%j),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     DO ist=1,staging_env%j-1
        iii(ist)=(ist-1)*staging_env%j+1 ! first el
        jjj(ist)=iii(ist)+staging_env%j-1 ! last el
        kkk(ist)=iii(ist)-1 ! prev el
     END DO
     kkk(1)=staging_env%p

     uf=f
     ! staging beads
     DO k=1,staging_env%nseg
        DO i=2,staging_env%j
           uf(i+iii(k),:)=uf(i+iii(k),:)&
                +REAL(i-1,dp)/REAL(i,dp)*uf(i+iii(k)-1,:)
        END DO
     END DO
     ! end point beads
     DO idim=1,SIZE(uf,2)
        DO k=1,staging_env%nseg
           sum_f=0._dp
           DO ij=2,staging_env%j-1
              sum_f=sum_f+uf((k-1)*staging_env%j+ij,idim)
           END DO
           uf(iii(k),idim)=uf(iii(k),idim)+&
                sum_f-const*(uf(jjj(k),idim)-uf(kkk(k),idim))
        END DO
     END DO
  END IF
END SUBROUTINE staging_f2uf
!***************************************************************************

!!****f* pint_types/pint_init_v *
!!
!!   NAME
!!     pint_init_v
!!
!!   FUNCTION
!!     initial velocity
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the pint env in which you should initialize the 
!!       velocity
!!     - kT: the temperature the particles should be initialized to
!!       (defaults to pint_env%kT)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_init_v(pint_env,kT,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), INTENT(in), OPTIONAL      :: kT
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: mykt

  failure=.FALSE.
  
  mykt=pint_env%kT
  IF (PRESENT(kT)) mykt=kT
  DO idim=1,SIZE(pint_env%uv,2)
     DO ib=1,SIZE(pint_env%uv,1)
        pint_env%uv(ib,idim) = next_random_number(rng_stream=pint_env%randomG,&
                                                  variance=kT/pint_env%mass_fict(ib,idim),&
                                                  error=error)
     END DO
  END DO
  DO idim=1,SIZE(pint_env%tv,3)
     DO ib=1,SIZE(pint_env%tv,2)
        DO inos=1,SIZE(pint_env%tv,1)
           pint_env%tv(inos,ib,idim) = next_random_number(rng_stream=pint_env%randomG,&
                                                          variance=kT/pint_env%Q(ib),&
                                                          error=error)
        END DO
     END DO
  END DO
END SUBROUTINE pint_init_v
!***************************************************************************

!!****f* pint_types/pint_calc_uf_h *
!!
!!   NAME
!!     pint_calc_uf_h
!!
!!   FUNCTION
!!     calculates the harmonic force in the u basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the path integral environment in which the harmonic
!!       forces should be calculated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION  pint_calc_uf_h(pint_env, error)RESULT(res)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_uf_h', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!  REAL(kind=dp), DIMENSION(:,:), INTENT(out), OPTIONAL, TARGET :: uf_h
!  REAL(kind=dp), DIMENSION(:,:), INTENT(in), OPTIONAL, TARGET :: ux

  failure=.FALSE.

  IF (.NOT. failure) THEN
     pint_env%e_pot_h=staging_calc_uf_h(pint_env%staging_env,&
          pint_env%mass_beads,pint_env%ux,pint_env%uf_h,error=error)
     res=pint_env%e_pot_h
     pint_env%uf_h=pint_env%uf_h/pint_env%mass_fict
  END IF
END FUNCTION pint_calc_uf_h
!***************************************************************************

!!****f* pint_types/staging_calc_uf_h *
!!
!!   NAME
!!     staging_calc_uf_h
!!
!!   FUNCTION
!!     calculates the harmonic force in the staging basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - staging_env: the staging environment
!!     - mass_beads: the masses of the beads
!!     - ux: the positions of the beads in the staging basis
!!     - uf_h: the harmonic forces (not accelerations)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION  staging_calc_uf_h(staging_env, mass_beads, ux, uf_h, error)RESULT(res)
    TYPE(staging_env_type), POINTER          :: staging_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in)                             :: mass_beads, ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out)                            :: uf_h
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'staging_calc_uf_h', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, isg, ist, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii, jjj, kkk
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: d, e_h, f

  failure=.FALSE.

  IF (.NOT. failure) THEN
     e_h=0._dp

     ALLOCATE(iii(staging_env%nseg),jjj(staging_env%nseg),&
          kkk(staging_env%nseg),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

     DO ist=1,staging_env%nseg
        iii(ist)=(ist-1)*staging_env%j+1 ! first el
        jjj(ist)=iii(ist)+staging_env%j  ! next fisrt (pbc)
        kkk(ist)=iii(ist)-staging_env%j  ! prev first el (pbc)
     END DO
     jjj(staging_env%nseg)=1
     kkk(1)=staging_env%p-staging_env%j

     DO idim=1,SIZE(mass_beads,2)
        DO ist=1,staging_env%nseg
           e_h=0.5*e_h+mass_beads(1,idim)*staging_env%w_j**2*&
                (ux(iii(ist),idim)-ux(jjj(ist),idim))**2
           uf_h(iii(ist),idim)=mass_beads(1,idim)*staging_env%w_j**2*(&
                2._dp*ux(iii(ist),idim)&
                -ux(jjj(ist),idim)&
                -ux(kkk(ist),idim)&
                )
           DO isg=2,staging_env%j ! use 3 as start?
              d=ux((ist-1)*staging_env%j+isg,idim)
              f=mass_beads((ist-1)*staging_env%j+isg,idim)*staging_env%w_j**2*d
              e_h=e_h+0.5_dp*f*d
              uf_h((ist-1)*staging_env%j+isg,idim)=f
           END DO
        END DO
     END DO
     res=e_h
  END IF
END FUNCTION staging_calc_uf_h
!***************************************************************************

!!****f* pint_types/pint_calc_f *
!!
!!   NAME
!!     pint_calc_f
!!
!!   FUNCTION
!!     calculates the force (and energy) in each bead, returns the sum
!!     of the potential energy
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: path integral environment on which you want calculate
!!       the forces
!!     - x: positions at which you want to evaluate the forces
!!     - f: the forces
!!     - e: potential energy on each bead
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_calc_f(pint_env,x,f,e,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: x
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: f
    REAL(kind=dp), DIMENSION(:), &
      INTENT(out), OPTIONAL                  :: e
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: res
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_f, my_x

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_x => pint_env%x
     IF (PRESENT(x)) my_x => x
     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           pint_env%replicas%r(idim,ib)=my_x(ib,idim)
        END DO
     END DO
     CALL rep_env_calc_e_f(pint_env%replicas,error=error)
     my_f => pint_env%f
     IF (PRESENT(f)) my_f => f
     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           my_f(ib,idim)=pint_env%replicas%f(idim,ib)
        END DO
     END DO
     IF (PRESENT(e)) THEN
        e=pint_env%replicas%f(SIZE(pint_env%replicas%f,2),:)
     ELSE
        pint_env%e_pot_bead = pint_env%replicas%f(SIZE(pint_env%replicas%f,2),:)
     END IF
     res=SUM(pint_env%replicas%f(:,SIZE(pint_env%replicas%f)))
!          /REAL(SIZE(pint_env%replicas%f),dp)
  END IF
END SUBROUTINE pint_calc_f
!***************************************************************************

!!****f* pint_types/pint_setup *
!!
!!   NAME
!!     pint_setup
!!
!!   FUNCTION
!!     prepares the things (forces) to perform a step of md
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: path integral env
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_setup(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: e_h

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pint_x2u(pint_env,error=error)
     e_h=pint_calc_uf_h(pint_env,error=error)
     CALL pint_calc_f(pint_env,error=error)
     CALL pint_f2uf(pint_env,error=error)
     DO idim=1,SIZE(pint_env%uf_h,2)
        DO ib=1,SIZE(pint_env%uf_h,1)
           pint_env%uf_h(ib,idim)=pint_env%uf_h(ib,idim)&
                +REAL(pint_env%nrespa,dp)*pint_env%uf(ib,idim)
        END DO
     END DO
     DO idim=1,SIZE(pint_env%uf_h,2)
        DO ib=1,SIZE(pint_env%uf_h,1)
           pint_env%tf(1,ib,idim)=(pint_env%mass_fict(ib,idim)*pint_env%uv(ib,idim)**2&
                -pint_env%kT)/pint_env%Q(ib)
        END DO
     END DO

     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           DO inos=1,pint_env%nnos-1
              pint_env%tf(inos+1,ib,idim)=pint_env%tv(inos,ib,idim)**2-&
                   pint_env%kT/pint_env%Q(ib)
           END DO
           DO inos=1,pint_env%nnos-1
              pint_env%tf(inos,ib,idim)=pint_env%tf(inos,ib,idim)&
                   -pint_env%tv(inos,ib,idim)*pint_env%tv(inos+1,ib,idim)
           END DO
        END DO
     END DO

  END IF
END SUBROUTINE pint_setup
!***************************************************************************

!!****f* pint_types/pint_calc_e_kin_beads_u *
!!
!!   NAME
!!     pint_calc_e_kin_beads_u
!!
!!   FUNCTION
!!     calculate the kinetic energy of the beads (in the u variables)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION pint_calc_e_kin_beads_u(pint_env,uv,error) RESULT(res)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: uv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_e_kin_beads_u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_uv

  failure=.FALSE.
  

  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  res=-1.0_dp
  IF (.NOT. failure) THEN
     IF (PRESENT(uv)) my_uv => uv
     res=0._dp
     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           res=res+pint_env%mass_fict(ib,idim)*my_uv(ib,idim)**2
        END DO
     END DO
     res=res*0.5
     IF (.not.PRESENT(uv)) pint_env%e_kin_beads=res
  END IF
END FUNCTION pint_calc_e_kin_beads_u
!***************************************************************************

!!****f* pint_types/pint_calc_nh_energy *
!!
!!   NAME
!!     pint_calc_nh_energy
!!
!!   FUNCTION
!!     calculates the energy (potential and kinetic) of the Nose-Hoover
!!     chain thermostats
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the path integral environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_calc_nh_energy(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_nh_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ekin, epot

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ekin=0._dp
     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           DO inos=1,pint_env%nnos
              ekin=ekin+pint_env%Q(ib)*pint_env%tv(inos,ib,idim)**2
           END DO
        END DO
     END DO
     pint_env%e_kin_t=0.5_dp*ekin
     epot=0._dp
     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           DO inos=1,pint_env%nnos
              epot=epot+pint_env%tx(inos,ib,idim)
           END DO
        END DO
     END DO
     pint_env%e_pot_t=pint_env%kT*epot
  END IF
END SUBROUTINE pint_calc_nh_energy
!***************************************************************************

!!****f* pint_types/pint_step *
!!
!!   NAME
!!     pint_step
!!
!!   FUNCTION
!!     does an md step (and nrespa harmonic evaluations)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: the path integral environment to evolve
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_step(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ib, idim, ierr, inos, iresp
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dti, dti2, dti22, e_h, rn, &
                                                tdti, tol
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: tmp
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env

  failure=.FALSE.

  IF (.NOT. failure) THEN
     rn=REAL(pint_env%nrespa,dp)
     dti=pint_env%dt/rn
     dti2=dti/2._dp
     tdti=2.*dti
     dti22=dti**2/2._dp
     DO iresp=1,pint_env%nrespa
        ! integrate bead positions
        IF (pint_env%nnos>1) THEN
           pint_env%ux=pint_env%ux+dti*pint_env%uv+dti22*(pint_env%uf_h-pint_env%uv*pint_env%tv(1,:,:))
        ELSE
           pint_env%ux=pint_env%ux+dti*pint_env%uv+dti22*pint_env%uf_h
        END IF
        ! integrate nh vars
        pint_env%tx=pint_env%tx+dti*pint_env%tv+dti22*pint_env%tf
        ! integrate v at half step
        IF (pint_env%nnos>1) THEN
           pint_env%uv_t=pint_env%uv+dti2*(pint_env%uf_h-pint_env%uv*pint_env%tv(1,:,:))
        ELSE
           pint_env%uv_t=pint_env%uv+dti2*pint_env%uf_h
        END IF

        tmp => pint_env%tv_t
        pint_env%tv_t => pint_env%tv
        pint_env%tv => tmp

        pint_env%tv=pint_env%tv_old+tdti*pint_env%tf
        pint_env%tv_old=pint_env%tv_t
        pint_env%tv_t=pint_env%tv_t+dti2*pint_env%tf

        ! calc forces at new pos
        e_h=pint_calc_uf_h(pint_env,error=error)
        IF (iresp==pint_env%nrespa) THEN
           CALL pint_u2x(pint_env,error=error)
           CALL pint_calc_f(pint_env,error=error)
           CALL pint_f2uf(pint_env,error=error)
           pint_env%uf_h=pint_env%uf_h+rn*pint_env%uf
        END IF

        ! add the new forces to v
        pint_env%uv_t=pint_env%uv_t+dti2*pint_env%uf_h
        ! iterate v to convergence
        IF (pint_env%nnos>0) THEN
           DO i=1,6
              tol=0._dp
              pint_env%uv_new=pint_env%uv_t/(1.+dti2*pint_env%tv(1,:,:))
              DO idim=1,pint_env%ndim
                 DO ib=1,pint_env%p
                    pint_env%tf(0,ib,idim)=(pint_env%mass_fict(ib,idim)*&
                         pint_env%uv_new(ib,idim)**2-pint_env%kT)/&
                         pint_env%Q(ib)
                 END DO
              END DO

              DO idim=1,pint_env%ndim
                 DO ib=1,pint_env%p
                    DO inos=1,pint_env%nnos-1
                       pint_env%tv_new(inos,:,:)=&
                            (pint_env%tv_t(inos,:,:)+dti2*pint_env%tf(inos,:,:))/&
                            (1._dp+dti2*pint_env%tv(inos+1,:,:))
                       pint_env%tf(inos+1,ib,idim)=&
                            (pint_env%tv_new(inos,ib,idim)**2-&
                            pint_env%kT/pint_env%Q(ib))
                       tol=MAX(tol,ABS(pint_env%tv(inos,ib,idim)&
                            -pint_env%tv_new(inos,ib,idim)))
                    END DO
                    pint_env%tv_new(pint_env%nnos,ib,idim)=&
                         pint_env%tv_t(pint_env%nnos,ib,idim)+&
                         dti2*pint_env%tf(pint_env%nnos,ib,idim)
                    tol=MAX(tol,ABS(pint_env%tv(pint_env%nnos,ib,idim)&
                         -pint_env%tv_new(pint_env%nnos,ib,idim)))
                    tol=MAX(tol,ABS(pint_env%uv(ib,idim)&
                         -pint_env%uv_new(ib,idim)))
                 END DO
              END DO

              pint_env%uv=pint_env%uv_new
              pint_env%tv=pint_env%tv_new
              IF (tol <= pint_env%v_tol) EXIT
           END DO
           CPPostcondition(tol>pint_env%v_tol,cp_warning_level,routineP,error,failure)
        ELSE
           pint_env%uv=pint_env%uv_t
        ENDIF
        DO inos=1,pint_env%nnos-1
           pint_env%tf(inos,:,:)=pint_env%tf(inos,:,:)&
                -pint_env%tv(inos,:,:)*pint_env%tv(inos+1,:,:)
        END DO
        CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
             f_env=f_env,new_error=new_error, failure=failure)
        logger => cp_error_get_logger(new_error)
        CALL cp_assert(logger%iter_info%iteration(2)==pint_env%iter,&
             cp_failure_level,cp_assertion_failed,routineP,&
             "md & force_eval lost sychro "//&
             CPSourceFileRef,&
             error,failure)
        CALL f_env_rm_defaults(f_env,new_error,ierr)
        pint_env%iter=pint_env%iter+1
        pint_env%t=pint_env%t+pint_env%dt     
     END DO
  END IF
END SUBROUTINE pint_step
!***************************************************************************

!!****f* pint_types/pint_write_centroids *
!!
!!   NAME
!!     pint_write_centroids
!!
!!   FUNCTION
!!     writes out the trajectory of the centroids (and their spread)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: path integral environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_write_centroids(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_centroids', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: ext, form, title
    INTEGER                                  :: i_form, iat, ib, idim, idir, &
                                                ierr, iter, unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: nb, ss, ss2
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: print_key

  failure=.FALSE.
  
  NULLIFY(f_env,cell,subsys)
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (pint_env%replicas%cart%mepos(2)==0) THEN
        print_key => section_vals_get_subs_vals(pint_env%input,&
             "MOTION%PINT%PRINT%CENTROIDS", error=error)
        CALL section_vals_val_get(print_key,"format",c_val=ext,i_val=i_form,&
             error=error)
        SELECT CASE(i_form)
        CASE (dump_dcd)
           form="UNFORMATTED"
        CASE default
           form="FORMATTED"
        END SELECT
        IF (BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
             basis_section=print_key, error=error),cp_p_file)) THEN
           unit_nr=cp_print_key_unit_nr(logger=logger, basis_section=print_key,&
                print_key_path="",extension="."//TRIM(ext), local=.FALSE., file_form=form, error=error)
           CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,f_env=f_env,&
                new_error=new_error, failure=failure)
           logger => cp_error_get_logger(new_error)
           IF (unit_nr>0) THEN
              iter=logger%iter_info%iteration(1)
              logger%iter_info%iteration(1)=0
              CALL force_env_get(force_env=f_env%force_env,&
                   cell=cell,&
                   subsys=subsys,&
                   error=new_error)
              CALL cp_subsys_get(subsys(1)%subsys,&
                   particles=particles,&
                   error=new_error)

              nb=REAL(pint_env%p,dp)
              idim=0
              DO iat=1,pint_env%ndim
                 DO idir=1,3
                    idim=idim+1
                    ss=0._dp
                    ss2=0._dp
                    DO ib=1,pint_env%p
                       ss=ss+pint_env%x(idim,ib)
                       ss2=ss2+pint_env%x(idim,ib)**2
                    END DO
                    particles%els(iat)%r(idir)=ss/nb
                    particles%els(iat)%v(idir)=SQRT(ss2/nb-(ss/nb)**2)
                 END DO
              END DO
              
              WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                   " i =",logger%iter_info%iteration(2),&
                   ", E =",SUM(pint_env%e_pot_bead)/REAL(pint_env%p,dp)
              CALL write_particle_coordinates(particles%els,cell=cell,&
                   iunit=unit_nr,FORMAT=i_form,content="POS_VEL",title=title)
              logger%iter_info%iteration(1)=iter
              CALL cp_print_key_finished_output(unit_nr,logger,&
                   print_key,"", local=.FALSE., error=error)
END IF
           CALL f_env_rm_defaults(f_env,new_error,ierr)
           CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
        END IF
     END IF
  END IF
END SUBROUTINE pint_write_centroids
!***************************************************************************

!!****f* pint_types/pint_write_ener *
!!
!!   NAME
!!     pint_write_ener
!!
!!   FUNCTION
!!     writes out the ener info of the centroids
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pint_env: path integral environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_write_ener(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iter, unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: cons, t, used_time
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (pint_env%replicas%cart%mepos(2)==0) THEN
        print_key => section_vals_get_subs_vals(pint_env%input,&
             "MOTION%PINT%PRINT%ENER_INFO", error=error)
        IF (BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
             basis_section=print_key, error=error),cp_p_file)) THEN
           unit_nr=cp_print_key_unit_nr(logger,print_key,"",&
                extension=".ener", local=.FALSE., error=error)
           IF (unit_nr>0) THEN
             ! please change also the corresponding format explaination above
             ! keep the constant of motion the true constant of motion !
             ! "Step Nr.    Time[fs]  Kin.[a.u.]     Temp[K]  Pot.[a.u.]    Cons Qty UsedTime[s]"
             !
              cons=(SUM(pint_env%e_pot_bead)+pint_env%e_pot_h+&
                   pint_env%e_kin_beads+pint_env%e_pot_t+pint_env%e_kin_t)
              t=cp_unit_from_cp2k(pint_env%t,"fs",error=error)
              used_time=0._dp ! to do
              WRITE (unit_nr,"(I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.9)")&
                   pint_env%iter,t,pint_env%e_kin_beads/REAL(pint_env%p,dp),&
                   SUM(pint_env%e_pot_bead)/REAL(pint_env%p,dp),&
                   cons,used_time
              CALL m_flush(unit_nr)
           END IF
           CALL cp_print_key_finished_output(unit_nr,logger,print_key,"",&
                local=.FALSE., error=error)
        END IF
     END IF
  END IF
END SUBROUTINE pint_write_ener
!***************************************************************************

!!****f* pint_types/pint_start [1.0] *
!!
!!   NAME
!!     pint_start
!!
!!   FUNCTION
!!     reads the restart, makes initial initialization to beads,...
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: parallel environment
!!     - input: the input to test
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pint_start(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_start', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, ii, inos, n_rep_val
    LOGICAL                                  :: centroid_speed, explicit, &
                                                failure, randomize_pos
    REAL(kind=dp)                            :: rv
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_vals
    TYPE(section_vals_type), POINTER         :: beads_section, coord_section, &
                                                input, nose_section, &
                                                velocity_section

  failure=.FALSE.
  
  NULLIFY(input, nose_section,coord_section,velocity_section)
  CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pint_init_v(pint_env,error=error)
     CALL section_vals_val_get(input,"MOTION%PINT%INIT%CENTROID_SPEED",&
          l_val=centroid_speed,error=error)
     IF (centroid_speed) THEN
        CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
        DO idim=1,pint_env%ndim
           rv = next_random_number(rng_stream=pint_env%randomG,&
                                   variance=pint_env%mass(idim)*pint_env%kT,&
                                   error=error)/pint_env%mass(idim)
           DO ib=1,pint_env%p
              pint_env%v(ib,idim)=pint_env%v(ib,idim)+rv
           END DO
        END DO
        CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
     END IF
     
     nose_section => section_vals_get_subs_vals(input,"MOTION%PINT%NOSE",&
          error=error)
     coord_section => section_vals_get_subs_vals(nose_section,"COORD",&
          error=error)
     CALL section_vals_get(coord_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
             n_rep_val=n_rep_val,error=error)
        IF (n_rep_val>0) THEN
           CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
           CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                i_rep_val=ib,r_vals=r_vals,error=error)
           CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim*pint_env%nnos,&
                cp_failure_level,cp_assertion_failed,"Invalid size of MOTION%PINT%NOSE%COORD "//&
                CPSourceFileRef,&
                routineP,error,failure)
           ii=0
           DO idim=1,pint_env%ndim
              DO ib=1,pint_env%p
                 DO inos=1,pint_env%nnos
                    ii=ii+1
                    pint_env%tx(inos,ib,idim)=r_vals(ii)
                 END DO
              END DO
           END DO
        END IF
     END IF
     velocity_section => section_vals_get_subs_vals(nose_section,"VELOCITY",&
          error=error)
     CALL section_vals_get(velocity_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(velocity_section,"_DEFAULT_KEYWORD_",&
             n_rep_val=n_rep_val,error=error)
        IF (n_rep_val>0) THEN
           CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
           CALL section_vals_val_get(velocity_section,"_DEFAULT_KEYWORD_",&
                i_rep_val=ib,r_vals=r_vals,error=error)
           CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim*pint_env%nnos,&
                cp_failure_level,cp_assertion_failed,"Invalid size of MOTION%PINT%NOSE%VELOCITY "//&
                CPSourceFileRef,&
                routineP,error,failure)
           ii=0
           DO idim=1,pint_env%ndim
              DO ib=1,pint_env%p
                 DO inos=1,pint_env%nnos
                    ii=ii+1
                    pint_env%tv(inos,ib,idim)=r_vals(ii)
                 END DO
              END DO
           END DO
        END IF
     END IF

     beads_section => section_vals_get_subs_vals(input,"MOTION%PINT%BEADS",&
          error=error)
     coord_section => section_vals_get_subs_vals(beads_section,"COORD",&
          error=error)
     CALL section_vals_get(coord_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
             n_rep_val=n_rep_val,error=error)
        IF (n_rep_val>0) THEN
           CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
           CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                i_rep_val=ib,r_vals=r_vals,error=error)
           CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim,cp_failure_level,&
                cp_assertion_failed,"Invalid size of MOTION%PINT%BEADS%COORD "//&
                CPSourceFileRef,&
                routineP,error,failure)
           ii=0
           DO idim=1,pint_env%ndim
              DO ib=1,pint_env%p
                    ii=ii+1
                    pint_env%x(ib,idim)=r_vals(ii)
              END DO
           END DO
        ELSE
           
        END IF
     END IF
     velocity_section => section_vals_get_subs_vals(beads_section,"VELOCITY",&
          error=error)
     CALL section_vals_get(velocity_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(velocity_section,"_DEFAULT_KEYWORD_",&
             n_rep_val=n_rep_val,error=error)
        IF (n_rep_val>0) THEN
           CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
           CALL section_vals_val_get(velocity_section,"_DEFAULT_KEYWORD_",&
                i_rep_val=ib,r_vals=r_vals,error=error)
           CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim,cp_failure_level,&
                cp_assertion_failed,"Invalid size of MOTION%PINT%BEAD%VELOCITY "//&
                CPSourceFileRef,&
                routineP,error,failure)
           ii=0
           DO idim=1,pint_env%ndim
              DO ib=1,pint_env%p
                 ii=ii+1
                 pint_env%v(ib,idim)=r_vals(ii)
              END DO
           END DO
        END IF
     END IF
     
     CALL section_vals_val_get(input,"MOTION%PINT%INIT%RANDOMIZE_POS",&
          l_val=randomize_pos,error=error)
     IF (randomize_pos) THEN
        DO idim=1,pint_env%ndim
           DO ib=1,pint_env%p
              pint_env%x(ib,idim) = pint_env%x(ib,idim) +&
                                    next_random_number(rng_stream=pint_env%randomG,&
                                                       variance=pint_env%beta/&
                                                                &SQRT(12.0_dp*pint_env%mass(idim)),&
                                    error=error)
           END DO
        END DO
     END IF
  END IF

  CALL pint_setup(pint_env,error=error)
  CALL pint_write_restart(pint_env,error=error)

END SUBROUTINE pint_start
!***************************************************************************

!!****f* pint_types/pint_write_restart [1.0] *
!!
!!   NAME
!!     pint_write_restart
!!
!!   FUNCTION
!!     writes a restart file for pint
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pint_write_restart(pint_env,force,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    LOGICAL, INTENT(in), OPTIONAL            :: force
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_write_restart', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     ! to do
  END IF
END SUBROUTINE pint_write_restart
!***************************************************************************

!!****f* pint_types/pint_run [1.0] *
!!
!!   NAME
!!     pint_run
!!
!!   FUNCTION
!!     runs the dynamic
!!
!!   NOTES
!!    everything should be read for an md step
!!
!!   INPUTS
!!     - pint_env: the pint env to run
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pint_run(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, num_steps
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: scal, t_tol

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL section_vals_val_get(pint_env%input,"MOTION%PINT%NUM_STEPS",&
          i_val=num_steps,error=error)
     CALL section_vals_val_get(pint_env%input,"MOTION%PINT%T_TOL",&
          i_val=num_steps,error=error)
     DO i=1,num_steps
        CALL pint_write_ener(pint_env,error=error)
        CALL pint_write_centroids(pint_env,error=error)
        CALL pint_write_restart(pint_env,error=error)
        IF (t_tol>0._dp) THEN
           IF (ABS(2._dp*pint_env%e_kin_beads/(pint_env%p*pint_env%ndim)&
                -pint_env%kT)>t_tol) THEN
              scal=SQRT(pint_env%kT*(pint_env%p*pint_env%ndim)/(2.*pint_env%e_kin_beads))
              pint_env%uv=scal*pint_env%uv
              CALL pint_setup(pint_env,error=error)
           END IF
        END IF
        CALL pint_step(pint_env,error=error)
     END DO
     ! set last flag!!!
     CALL pint_write_restart(pint_env,error=error)
  END IF
END SUBROUTINE pint_run
!***************************************************************************

END MODULE pint_types
