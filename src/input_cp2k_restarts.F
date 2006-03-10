!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/input_cp2k_restarts [1.0] *
!!
!!   NAME
!!     input_cp2k_restarts
!!
!!   FUNCTION
!!     Set of routines to dump the restart file of CP2K
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     01.2006 [created] Teodoro Laino     
!! 
!!
!!   SOURCE
!******************************************************************************

MODULE input_cp2k_restarts
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_linked_list_val,              ONLY: cp_sll_val_create,&
                                             cp_sll_val_get_length,&
                                             cp_sll_val_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: use_unit_angstrom,&
                                             use_unit_bohr,&
                                             use_unit_scaled_angstrom,&
                                             use_unit_scaled_bohr
  USE input_section_types,             ONLY: section_get_keyword_index,&
                                             section_type,&
                                             section_vals_add_values,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_val_unset
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_sum
  USE metadynamics_types,              ONLY: meta_env_type
  USE parallel_rng_types,              ONLY: dump_rng_stream,&
                                             rng_record_length
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: bohr
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: ascii_to_string,&
                                             string_to_ascii
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_restarts'

  PUBLIC :: update_input            
  
!!*****
!******************************************************************************

CONTAINS

! *****************************************************************************

!!****f* input_cp2k_restarts/update_input [1.0] *
!!
!!   NAME
!!     update_input
!!
!!   FUNCTION
!!     Updates the whole input file for the restart
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE update_input(md_env,force_env,kg_env,meta_env,globenv,force_env_section,error)
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_input', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, lcond
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(force_env_type), POINTER            :: my_force_env
    TYPE(kg_environment_type), POINTER       :: my_kg_env
    TYPE(meta_env_type), POINTER             :: my_meta_env
    TYPE(section_vals_type), POINTER         :: motion_section

    failure = .FALSE.
    NULLIFY(logger,motion_section,my_force_env,my_kg_env,my_meta_env,force_env_p)
    logger => cp_error_get_logger(error)
    !
    ! Can handle md_env or force_env
    !
    lcond = PRESENT(md_env).OR.PRESENT(force_env).OR.PRESENT(kg_env)
    CPPrecondition(lcond,cp_failure_level,routineP,error,failure)
    IF (PRESENT(md_env)) THEN
       CALL get_md_env ( md_env=md_env,&
                         force_env_p=force_env_p)

       CPPrecondition(ASSOCIATED(force_env_p),cp_failure_level,routineP,error,failure)
       my_force_env => force_env_p(1)%force_env       
    ELSE IF (PRESENT(force_env)) THEN 
       my_force_env => force_env
    END IF
    !
    ! Let's set up kg_env
    !
    IF (PRESENT(kg_env)) THEN
       my_kg_env => kg_env
    ELSEIF (ASSOCIATED(my_force_env)) THEN 
       my_kg_env => my_force_env%kg_env
    END IF
    !
    ! meta_env
    !
    IF (PRESENT(meta_env)) THEN
       my_meta_env => meta_env
    ELSEIF (ASSOCIATED(my_force_env)) THEN 
       my_meta_env => my_force_env%meta_env
    END IF    
    !
    ! The real restart setting...
    !
    motion_section => section_vals_get_subs_vals(globenv%input_file,"MOTION",error=error)
    CALL update_motion(motion_section,md_env,my_force_env,my_kg_env,my_meta_env,&
         logger,error)
    !
    ! Here we have to handle multiple force_env, in order to map them to the 
    ! proper force_eval section.. To fix for the future..
    !
    IF (ASSOCIATED(my_force_env)) &
         CALL update_force_eval(force_env_section,my_force_env,logger,error)

  END SUBROUTINE update_input

!!****f* input_cp2k_restarts/update_motion [1.0] *
!!
!!   NAME
!!     update_motion
!!
!!   FUNCTION
!!     Updates the motion section of the input file
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE update_motion(motion_section,md_env,force_env,kg_env,meta_env,logger,error)
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_motion', &
      routineP = moduleN//':'//routineN

    INTEGER :: counter, handle, i, iatom, icoef, iel, ikind, ioff, iproc, &
      irep, j, my_accepted_md, my_itimes, my_md_approach_nr, n_el, natom, &
      nhc_len, nkind, nsgf, nsgf_max, num_nhc, numneed, offset, stat, &
      tot_nhcneed
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cindex
    INTEGER, POINTER                         :: accepted_md, INDEX(:), &
                                                itimes, md_approach_nr, &
                                                nhc_per_proc( : )
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_t
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbuffer, vbuffer
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer, eta, fnhc, mnhc, &
                                                veta, work
    REAL(KIND=dp), POINTER                   :: t
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(lnhc_parameters_type), POINTER      :: nhc_baro( : ), nhc_coef( : ), &
                                                nhc_part( : )
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    CALL timeset(routineN,"I",'',handle)
    NULLIFY (para_env,logger,nhc_part,nhc_baro,npt,para_env)
    NULLIFY (atomic_kind_set, work_section, dyn_coeff_set, md_approach_nr, accepted_md)
    NULLIFY (atomic_kind, coeffs, aux_basis_set, nhc_coef, t, itimes)
    IF (PRESENT(md_env)) THEN
       CALL get_md_env ( md_env=md_env,&
                         nhc_part=nhc_part,&
                         nhc_baro=nhc_baro,&
                         nhc_coef=nhc_coef,&
                         npt=npt,&
                         t=t,&
                         itimes=itimes,&
                         accepted_md=accepted_md,&
                         md_approach_nr=md_approach_nr,&
                         para_env=para_env)
    ELSE
       IF (ASSOCIATED(force_env)) THEN
          para_env => force_env%para_env
       ELSEIF (ASSOCIATED(kg_env)) THEN
          para_env => kg_env%para_env
       ELSE
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    
    IF (ASSOCIATED(itimes)) THEN
       IF (itimes>=0) THEN
          my_itimes = itimes
          CALL section_vals_val_set(motion_section,"MD%STEP_START_VAL",i_val=my_itimes,error=error)
          CPAssert(ASSOCIATED(t),cp_failure_level,routineP,error,failure)
          my_t = t
          CALL section_vals_val_set(motion_section,"MD%TIME_START_VAL",r_val=my_t,error=error)
       END IF
    END IF
    IF (ASSOCIATED(md_approach_nr)) THEN
       IF (md_approach_nr>=0) THEN
          my_md_approach_nr = md_approach_nr
          CALL section_vals_val_set(motion_section,"MD%HMC%MD_APPROACH_START_VAL",i_val=my_md_approach_nr,&
               error=error)
          CPAssert(ASSOCIATED(accepted_md),cp_failure_level,routineP,error,failure)
          my_accepted_md = accepted_md
          CALL section_vals_val_set(motion_section,"MD%HMC%ACCEPTED_MD_START_VAL",i_val=my_accepted_md,&
               error=error)
       END IF
    END IF
    !
    ! PARTICLE THERMOSTAT
    !
    IF (ASSOCIATED(nhc_part)) THEN
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! write out nhc_part variable in a consistent way for restarts
       !-----------------------------------------------------------------------------
       nhc_len = SIZE(nhc_part(1)%nvt,1)
       num_nhc = SIZE(nhc_part(1)%nvt,2)
       numneed = num_nhc
       ALLOCATE ( nhc_per_proc( para_env % num_pe ), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       nhc_per_proc ( : ) = 0

       CALL mp_allgather(numneed,nhc_per_proc,para_env%group)

       IF (nhc_part ( 1 ) % dis_type=='INTER_REP') THEN
          tot_nhcneed = nhc_per_proc(1)
       ELSE
          tot_nhcneed = SUM ( nhc_per_proc )
       ENDIF

       NULLIFY(work,index)

       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! nhc_part%eta  
       !-----------------------------------------------------------------------------
       ALLOCATE(eta(tot_nhcneed*nhc_len), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iproc=1,para_env%num_pe
          CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
          CALL reallocate(index,1,nhc_per_proc(iproc))
          IF(para_env%mepos == (iproc-1)) THEN
             INDEX(:)=0
             counter=0
             DO i=1,nhc_len
                DO j=1,num_nhc
                   counter=counter+1
                   work(counter)=nhc_part(1)%nvt(i,j)%eta
                   INDEX(j) = nhc_part ( 1 ) % INDEX ( j )
                END DO
             END DO
          ELSE
             work(:) = 0.0_dp
          END IF
          CALL mp_bcast(work,iproc-1,para_env%group)
          CALL mp_bcast(index,iproc-1,para_env%group)
          counter=0
          DO i=1,nhc_len
             DO j=1,nhc_per_proc(iproc)
                counter = counter+1
                eta((INDEX(j)-1)*nhc_len+i) = work(counter)
             END DO
          END DO
       END DO


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! nhc_part % veta   
       !-----------------------------------------------------------------------------
       ALLOCATE(veta(tot_nhcneed*nhc_len), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iproc=1,para_env%num_pe
          CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
          CALL reallocate(index,1,nhc_per_proc(iproc))
          IF(para_env%mepos == (iproc-1)) THEN
             INDEX(:)=0
             counter=0
             DO i=1,nhc_len
                DO j=1,num_nhc
                   counter=counter+1
                   work(counter)=nhc_part(1)%nvt(i,j)%v
                   INDEX(j) = nhc_part ( 1 ) % index ( j )
                END DO
             END DO
          ELSE
             work(:) = 0.0_dp
          END IF
          CALL mp_bcast(work,iproc-1,para_env%group)
          CALL mp_bcast(index,iproc-1,para_env%group)
          counter=0
          DO i=1,nhc_len
             DO j=1,nhc_per_proc(iproc)
                counter = counter+1
                veta((INDEX(j)-1)*nhc_len+i) = work(counter)
             END DO
          END DO
       END DO

       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! nhc_part % force
       !-----------------------------------------------------------------------------
       ALLOCATE(fnhc(tot_nhcneed*nhc_len), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iproc=1,para_env%num_pe
          CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
          CALL reallocate(index,1,nhc_per_proc(iproc))
          IF(para_env%mepos == (iproc-1)) THEN
             INDEX(:)=0
             counter=0
             DO i=1,nhc_len
                DO j=1,num_nhc
                   counter=counter+1
                   work(counter)=nhc_part(1)%nvt(i,j)%f
                   INDEX(j) = nhc_part ( 1 ) % index ( j )
                END DO
             END DO
          ELSE
             work(:) = 0.0_dp
          END IF
          CALL mp_bcast(work,iproc-1,para_env%group)
          CALL mp_bcast(index,iproc-1,para_env%group)
          counter=0
          DO i=1,nhc_len
             DO j=1,nhc_per_proc(iproc)
                counter = counter+1
                fnhc((INDEX(j)-1)*nhc_len+i) = work(counter)
             END DO
          END DO
       END DO


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! nhc_part % mass
       !-----------------------------------------------------------------------------
       ALLOCATE(mnhc(tot_nhcneed*nhc_len), STAT = stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)    
       DO iproc=1,para_env%num_pe
          CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
          CALL reallocate(index,1,nhc_per_proc(iproc))
          IF(para_env%mepos == (iproc-1)) THEN
             INDEX(:)=0
             counter=0
             DO i=1,nhc_len
                DO j=1,num_nhc
                   counter=counter+1
                   work(counter)=nhc_part(1)%nvt(i,j)%mass
                   INDEX(j) = nhc_part ( 1 ) % index ( j )
                END DO
             END DO
          ELSE
             work(:) = 0.0_dp
          END IF
          CALL mp_bcast(work,iproc-1,para_env%group)
          CALL mp_bcast(index,iproc-1,para_env%group)
          counter=0
          DO i=1,nhc_len
             DO j=1,nhc_per_proc(iproc)
                counter = counter+1
                mnhc((INDEX(j)-1)*nhc_len+i) = work(counter)
             END DO
          END DO
       END DO

       DEALLOCATE ( work, STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       DEALLOCATE ( index, STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       DEALLOCATE ( nhc_per_proc, STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 

       IF ( nhc_part ( 1 ) % dis_type=='INTER_REP' ) THEN
          eta  = 0.0_dp
          veta = 0.0_dp
          counter = 0
          DO i=1,SIZE(nhc_part(1)%nvt,1)
             DO j=1,SIZE(nhc_part(1)%nvt,2)
                counter = counter + 1
                eta (counter)= nhc_part(1)%nvt(i,j)%eta
                veta(counter)= nhc_part(1)%nvt(i,j)%v
                fnhc(counter)= nhc_part(1)%nvt(i,j)%f
                mnhc(counter)= nhc_part(1)%nvt(i,j)%mass
             END DO
          END DO
       END IF

       work_section => section_vals_get_subs_vals(motion_section,"MD%NOSE",error=error)
       CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)
    END IF
    !
    ! BAROSTAT
    !
    IF (ASSOCIATED(nhc_baro)) THEN
       nhc_len     = SIZE(nhc_baro(1)%nvt,1)
       tot_nhcneed = SIZE(nhc_baro(1)%nvt,2)
       ALLOCATE(eta(tot_nhcneed*nhc_len), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(veta(tot_nhcneed*nhc_len), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(fnhc(tot_nhcneed*nhc_len), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mnhc(tot_nhcneed*nhc_len), STAT = stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       
       counter = 0
       DO i=1,SIZE(nhc_baro(1)%nvt,1)
          DO j=1,SIZE(nhc_baro(1)%nvt,2)
             counter = counter + 1
             eta (counter)= nhc_baro(1)%nvt(i,j)%eta
             veta(counter)= nhc_baro(1)%nvt(i,j)%v
             fnhc(counter)= nhc_baro(1)%nvt(i,j)%f
             mnhc(counter)= nhc_baro(1)%nvt(i,j)%mass
          END DO
       END DO
       
       work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%NOSE",error=error)
       CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)
       
    END IF
    IF (ASSOCIATED(npt)) THEN
       ALLOCATE(veta(SIZE(npt,1)*SIZE(npt,2)), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mnhc(SIZE(npt,1)*SIZE(npt,2)), STAT = stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       
       counter = 0
       DO i=1,SIZE(npt,1)
          DO j=1,SIZE(npt,2)
             counter = counter + 1
             veta(counter) = npt(i,j)%v
             mnhc(counter) = npt(i,j)%mass
          END DO
       END DO
       work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT",error=error)
       CALL set_template_restart(work_section, veta=veta, mnhc=mnhc, error=error)
    END IF
    !
    ! COEFFICIENTS 
    !
    IF (ASSOCIATED(kg_env)) THEN
       CALL get_kg_env(kg_env=kg_env,&
            atomic_kind_set=atomic_kind_set,&
            dyn_coeff_set=dyn_coeff_set)
       
       IF (ASSOCIATED(dyn_coeff_set)) THEN
          
          nsgf_max=0
          n_el=0
          DO ikind=1,SIZE(atomic_kind_set)
             atomic_kind=> atomic_kind_set(ikind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  aux_basis_set=aux_basis_set,&
                  natom=natom)
             IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
             CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                  nsgf=nsgf)
             nsgf_max=MAX(nsgf_max,nsgf)
             n_el = n_el + natom
          END DO
          
          ALLOCATE (pbuffer(n_el,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (vbuffer(n_el,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (cindex(n_el),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          
          pbuffer(:,:)=0.0_dp
          vbuffer(:,:)=0.0_dp
          
          cindex(:)=0
          nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
          offset = 0
          DO ikind=1, nkind
             coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
             atomic_kind=> atomic_kind_set(ikind)
             IF(.NOT.ASSOCIATED(coeffs))THEN
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     natom=natom)
                offset = offset + natom
             ELSE
                DO iel= 1,coeffs%n_els
                   iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
                   ioff=iatom - offset
                   cindex(ioff)=iatom
                   DO icoef=1,coeffs%ncoef_atom
                      pbuffer(ioff,icoef)=coeffs%pos(iel,icoef)
                      vbuffer(ioff,icoef)=coeffs%vel(iel,icoef)
                   END DO
                END DO
             END IF
          END DO
          CALL mp_sum(pbuffer,para_env%group)
          CALL mp_sum(vbuffer,para_env%group)
          CALL mp_sum(cindex, para_env%group)
          
          ALLOCATE(INDEX(SIZE(cindex)),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          index = cindex
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%RESTART_INDEX",error=error)
          CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",i_vals_ptr=index,error=error)
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%COORD",error=error)
          DO irep = 1, SIZE(index)
             ALLOCATE(buffer(SIZE(pbuffer,2)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = pbuffer(irep,:)
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_vals_ptr=buffer,&
                  i_rep_val=irep,error=error)    
          END DO
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%VELOCITY",error=error)
          DO irep = 1, SIZE(index)
             ALLOCATE(buffer(SIZE(vbuffer,2)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = vbuffer(irep,:)       
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_vals_ptr=buffer,&
                  i_rep_val=irep,error=error)    
          END DO
          DEALLOCATE (pbuffer,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          DEALLOCATE (vbuffer,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          DEALLOCATE (cindex,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(nhc_coef)) THEN
          !
          ! NOSE OF KG_COEFF
          !
          nhc_len     = SIZE(nhc_coef(1)%nvt,1)
          tot_nhcneed = SIZE(nhc_coef(1)%nvt,2)
          ALLOCATE(eta(tot_nhcneed*nhc_len), STAT = stat )
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(veta(tot_nhcneed*nhc_len), STAT = stat )
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fnhc(tot_nhcneed*nhc_len), STAT = stat )
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mnhc(tot_nhcneed*nhc_len), STAT = stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
          
          counter = 0
          DO i=1,SIZE(nhc_coef(1)%nvt,1)
             DO j=1,SIZE(nhc_coef(1)%nvt,2)
                counter = counter + 1
                eta (counter)= nhc_coef(1)%nvt(i,j)%eta
                veta(counter)= nhc_coef(1)%nvt(i,j)%v
                fnhc(counter)= nhc_coef(1)%nvt(i,j)%f
                mnhc(counter)= nhc_coef(1)%nvt(i,j)%mass
             END DO
          END DO
          
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%NOSE",error=error)
          CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)    
       END IF
    END IF

    IF (ASSOCIATED(meta_env)) THEN
       CALL section_vals_val_set(meta_env%metadyn_section,"STEP_START_VAL",&
            i_val=meta_env%n_steps,error=error)
       CALL section_vals_val_set(meta_env%metadyn_section,"NHILLS_START_VAL",&
            i_val=meta_env%hills_env%n_hills,error=error)
       IF (meta_env%do_hills.AND.meta_env%hills_env%n_hills/=0) THEN
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_POS",error=error)
          CALL meta_hills_val_set_ss(work_section, meta_env, error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_AMPL",error=error)
          CALL meta_hills_val_set_ds(work_section, meta_env, error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_HEIGHT",error=error)
          CALL meta_hills_val_set_ww(work_section, meta_env, error=error)
       END IF
       IF (meta_env%extended_lagrange) THEN
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS0",error=error)
          DO irep = 1, meta_env%n_colvar
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_val=meta_env%metavar(irep)%ss0,&
                  i_rep_val=irep,error=error) 
          END DO
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_VVP",error=error)
          DO irep = 1, meta_env%n_colvar
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_val=meta_env%metavar(irep)%vvp,&
                  i_rep_val=irep,error=error) 
          END DO
       END IF
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE update_motion

!!****f* input_cp2k_restarts/set_template_restart [1.0] *
!!
!!   NAME
!!     set_template_restart
!!
!!   FUNCTION
!!     Set the nose structure like restart
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE set_template_restart(work_section, eta, veta, fnhc, mnhc, error)
    TYPE(section_vals_type), POINTER         :: work_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eta, veta, fnhc, mnhc
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_template_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: coord, force, mass, velocity

    failure = .FALSE.
    NULLIFY(coord, force, velocity, mass)
    IF (PRESENT(eta)) THEN
       IF (SIZE(eta)>0) THEN
          coord => section_vals_get_subs_vals(work_section,"COORD",error=error)
          CALL section_vals_val_set(coord,"_DEFAULT_KEYWORD_",r_vals_ptr=eta,error=error)
       ELSE
          DEALLOCATE(eta, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
    IF (PRESENT(veta)) THEN
       IF (SIZE(veta)>0) THEN
          velocity => section_vals_get_subs_vals(work_section,"VELOCITY",error=error)
          CALL section_vals_val_set(velocity,"_DEFAULT_KEYWORD_",r_vals_ptr=veta,error=error)
       ELSE
          DEALLOCATE(veta, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
    IF (PRESENT(fnhc)) THEN
       IF (SIZE(fnhc)>0) THEN
          force => section_vals_get_subs_vals(work_section,"FORCE",error=error)
          CALL section_vals_val_set(force,"_DEFAULT_KEYWORD_",r_vals_ptr=fnhc,error=error)
       ELSE
          DEALLOCATE(fnhc, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
    IF (PRESENT(mnhc)) THEN
       IF (SIZE(mnhc)>0) THEN
          mass => section_vals_get_subs_vals(work_section,"MASS",error=error)
          CALL section_vals_val_set(mass,"_DEFAULT_KEYWORD_",r_vals_ptr=mnhc,error=error)
       ELSE
          DEALLOCATE(mnhc, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
  END SUBROUTINE set_template_restart

!!****f* input_cp2k_restarts/update_force_eval [1.0] *
!!
!!   NAME
!!     update_force_eval
!!
!!   FUNCTION
!!     Updates the force_eval section of the input file
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE update_force_eval(force_env_section,force_env,logger,error)
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_force_eval', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: rng_section, subsys_section

    NULLIFY(rng_section, subsys_section)

    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    rng_section    => section_vals_get_subs_vals(subsys_section,"RNG_INIT",error=error)
    CALL update_rng(rng_section,force_env,error)
    CALL update_subsys(subsys_section,force_env,error)

  END SUBROUTINE update_force_eval

!!****f* input_cp2k_restarts/update_rng [1.0] *
!!
!!   NAME
!!     update_rng
!!
!!   FUNCTION
!!     Updates the rng section of the input file
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!*** **********************************************************************

  SUBROUTINE update_rng(rng_section,force_env,error)

    ! Write current status of the parallel random number generator (RNG)

    TYPE(section_vals_type), POINTER         :: rng_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_rng', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, isubsys, &
      nparticle, nparticle_kind, nparticle_local
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ascii
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles

    ! -------------------------------------------------------------------------

    CALL force_env_get(force_env, subsys=subsys, para_env=para_env, error=error)

    ! The structure for multiple subsys has to be still decided

    CPPostcondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                         atomic_kinds=atomic_kinds,&
                         local_particles=local_particles,&
                         particles=particles,&
                         error=error)

      IF (.NOT.ASSOCIATED(local_particles%local_particle_set)) CYCLE

      nparticle_kind = atomic_kinds%n_els
      nparticle = particles%n_els

      ALLOCATE (ascii(rng_record_length,nparticle),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ascii = 0

      DO iparticle=1,nparticle
        DO iparticle_kind=1,nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
            IF (iparticle == local_particles%list(iparticle_kind)%array(iparticle_local)) THEN
              CALL dump_rng_stream(rng_stream=local_particles%local_particle_set(iparticle_kind)%&
                                              &rng(iparticle_local)%stream,&
                                   rng_record=rng_record)
              CALL string_to_ascii(rng_record,ascii(:,iparticle))
            END IF
          END DO
        END DO
      END DO

      CALL mp_sum(ascii,para_env%group)

      CALL section_rng_val_set(rng_section=rng_section,&
                               particles=particles,&
                               ascii=ascii,&
                               error=error)

      DEALLOCATE (ascii,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO

  END SUBROUTINE update_rng

  ! ***************************************************************************

!!****f* input_cp2k_restarts/update_subsys [1.0] *
!!
!!   NAME
!!     update_subsys
!!
!!   FUNCTION
!!     Updates the subsys section of the input file
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE update_subsys(subsys_section,force_env,error)
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_subsys', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, isubsys, my_unit, &
                                                stat
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp)                            :: conv_factor
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    NULLIFY(work, work_section, particles, subsys, cell, cell_ref)
    CALL timeset(routineN,"I",'',handle)
    CALL force_env_get(force_env,&
                       subsys=subsys,&
                       cell=cell,&
                       cell_ref=cell_ref,&
                       use_ref_cell=use_ref_cell,&
                       error=error)
    !
    ! The structure for multiple subsys has to be still decided..
    !
    CPPostcondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
    DO isubsys=1,SIZE(subsys)
       CALL cp_subsys_get(subsys(isubsys)%subsys,&
            particles=particles, error=error)       
       work_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
       CALL section_vals_val_get(subsys_section,"CELL%UNITS",i_val=my_unit,error=error)
       SELECT CASE(my_unit)
       CASE (use_unit_angstrom,use_unit_scaled_angstrom)
          conv_factor = bohr
       CASE (use_unit_bohr,use_unit_scaled_bohr)
          conv_factor = 1.0_dp
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       !    Coordinates and Velocities
       work_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
       CALL section_coord_val_set(work_section, particles, conv_factor, error)
       work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
       CALL section_velocity_val_set(work_section, particles, conv_factor=1.0_dp, error=error)
       !
       ! CELL
       ! 
       work_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !     CELL VECTORS - A
       work(1:3) = cell%hmat(1,1:3)/conv_factor
       CALL section_vals_val_set(work_section,"A",r_vals_ptr=work,error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !     CELL VECTORS - B
       work(1:3) = cell%hmat(2,1:3)/conv_factor
       CALL section_vals_val_set(work_section,"B",r_vals_ptr=work,error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !     CELL VECTORS - C
       work(1:3) = cell%hmat(3,1:3)/conv_factor
       CALL section_vals_val_set(work_section,"C",r_vals_ptr=work,error=error)
       CALL section_vals_val_unset(work_section,"ABC",error=error)
       !
       ! CELL_REF
       !
       IF (use_ref_cell.AND.ASSOCIATED(cell_ref)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"CELL%CELL_REF",error=error)
          CALL section_vals_val_set(work_section,"UNITS",i_val=my_unit,error=error)
          ALLOCATE(work(3),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !     CELL VECTORS - A
          work(1:3) = cell_ref%hmat(1,1:3)/conv_factor
          CALL section_vals_val_set(work_section,"A",r_vals_ptr=work,error=error)
          ALLOCATE(work(3),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !     CELL VECTORS - B
          work(1:3) = cell_ref%hmat(2,1:3)/conv_factor
          CALL section_vals_val_set(work_section,"B",r_vals_ptr=work,error=error)
          ALLOCATE(work(3),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !     CELL VECTORS - C
          work(1:3) = cell_ref%hmat(3,1:3)/conv_factor
          CALL section_vals_val_set(work_section,"C",r_vals_ptr=work,error=error)
       END IF
    END DO
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE update_subsys


!!****f* input_cp2k_restarts/section_coord_val_set [1.0] *
!!
!!   NAME
!!     section_coord_val_set
!!
!!   FUNCTION
!!     routine to dump coordinates.. fast implementation
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE section_coord_val_set(coord_section, particles, conv_factor, error)

    TYPE(section_vals_type), POINTER         :: coord_section
    TYPE(particle_list_type), POINTER        :: particles
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_coord_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: line
    CHARACTER(LEN=default_string_length)     :: my_tag, name, rcval
    INTEGER                                  :: ik, irk, k, Nlist
    LOGICAL                                  :: failure, ldum
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(coord_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(coord_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => coord_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(coord_section%values,2)==1) EXIT
       CALL section_vals_add_values(coord_section,error=error)
    END DO
    vals  => coord_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,particles%n_els
       CALL get_atomic_kind(particles%els(irk)%atomic_kind,name=name)
       my_tag = particles%els(irk)%molname
       ldum = qmmm_ff_precond_only_qm(name)
       ldum = qmmm_ff_precond_only_qm(my_tag)
       line = TRIM(name)
       DO k = 1, 3
          WRITE(rcval,"(ES25.16)") particles%els(irk)%r(k)/conv_factor
          line = TRIM(line)//" "//TRIM(rcval)
       END DO
       line = TRIM(line)//" "//TRIM(my_tag)
       CALL val_create(my_val,lc_val=line,error=error)
       
       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos         
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    coord_section%values(ik,1)%list => vals
  END SUBROUTINE section_coord_val_set

!!****f* input_cp2k_restarts/section_velocity_val_set [1.0] *
!!
!!   NAME
!!     section_velocity_val_set
!!
!!   FUNCTION
!!     routine to dump velocities.. fast implementation
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE section_velocity_val_set(velocity_section, particles, conv_factor, error)

    TYPE(section_vals_type), POINTER         :: velocity_section
    TYPE(particle_list_type), POINTER        :: particles
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_velocity_val_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(velocity_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(velocity_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => velocity_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(velocity_section%values,2)==1) EXIT
       CALL section_vals_add_values(velocity_section,error=error)
    END DO
    vals  => velocity_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,particles%n_els
       ALLOCATE(vel(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       vel = particles%els(irk)%v(1:3)/conv_factor       
       CALL val_create(my_val,r_vals_ptr=vel,error=error)
       
       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos         
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    velocity_section%values(ik,1)%list => vals
  END SUBROUTINE section_velocity_val_set

!!****f* input_cp2k_restarts/section_rng_val_set [1.0] *
!!
!!   NAME
!!     section_rng_val_set
!!
!!   FUNCTION
!!     routine to dump rngs.. fast implementation
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [teo]
!!     - string dump (again) instead of integer ASCII code (07.03.06,MK)
!!
!!*** **********************************************************************

  SUBROUTINE section_rng_val_set(rng_section, particles, ascii, error)

    TYPE(section_vals_type), POINTER         :: rng_section
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, DIMENSION(:, :)                 :: ascii
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_rng_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: ik, irk, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(rng_section),cp_failure_level,routineP,error,failure)
    CPPrecondition((rng_section%ref_count > 0),cp_failure_level,routineP,error,failure)

    NULLIFY (my_val,old_val,section,vals)

    section => rng_section%section

    ik = section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)

    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
                   "section "//TRIM(section%name)//" does not contain keyword "//&
                   "_DEFAULT_KEYWORD_",error,failure)

    DO
      IF (SIZE(rng_section%values,2)==1) EXIT
      CALL section_vals_add_values(rng_section,error=error)
    END DO

    vals => rng_section%values(ik,1)%list
    Nlist = 0

    IF (ASSOCIATED(vals)) THEN
      Nlist = cp_sll_val_get_length(vals,error)
    END IF

    DO irk=1,particles%n_els

      CALL ascii_to_string(ascii(:,irk),rng_record)
      CALL val_create(val=my_val,lc_val=rng_record,error=error)

      IF (Nlist /= 0) THEN
        IF (irk == 1) THEN
          new_pos => vals
        ELSE
          new_pos => new_pos%rest
        END IF
        old_val => new_pos%first_el
        CALL val_release(old_val,error=error)
        new_pos%first_el => my_val
      ELSE
        IF (irk == 1) THEN
          NULLIFY (new_pos)
          CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
          vals => new_pos
        ELSE
          NULLIFY (new_pos%rest)
          CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
          new_pos => new_pos%rest
        END IF
      END IF

      NULLIFY (my_val)

    END DO

    rng_section%values(ik,1)%list => vals

  END SUBROUTINE section_rng_val_set

!!****f* input_cp2k_restarts/meta_hills_val_set_ss [1.0] *
!!
!!   NAME
!!     meta_hills_val_set_ss
!!
!!   FUNCTION
!!     routine to dump hills information during metadynamics run
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE meta_hills_val_set_ss(ss_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: ss_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_hills_val_set_ss', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ss_val
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(ss_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ss_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => ss_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(ss_section%values,2)==1) EXIT
       CALL section_vals_add_values(ss_section,error=error)
    END DO
    vals  => ss_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = SIZE(meta_env%hills_env%ss_history,1)
    DO irk=1,meta_env%hills_env%n_hills
       ALLOCATE(ss_val(lsize),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       ss_val = meta_env%hills_env%ss_history(:,irk) 
       CALL val_create(my_val,r_vals_ptr=ss_val,error=error)
       
       IF (irk <= Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos         
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    ss_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_ss

!!****f* input_cp2k_restarts/meta_hills_val_set_ds [1.0] *
!!
!!   NAME
!!     meta_hills_val_set_ds
!!
!!   FUNCTION
!!     routine to dump hills information during metadynamics run
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE meta_hills_val_set_ds(ds_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: ds_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_hills_val_set_ds', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ds_val
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(ds_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ds_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => ds_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(ds_section%values,2)==1) EXIT
       CALL section_vals_add_values(ds_section,error=error)
    END DO
    vals  => ds_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = SIZE(meta_env%hills_env%delta_s_history,1)
    DO irk=1,meta_env%hills_env%n_hills
       ALLOCATE(ds_val(lsize),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       ds_val = meta_env%hills_env%delta_s_history(:,irk) 
       CALL val_create(my_val,r_vals_ptr=ds_val,error=error)
       
       IF (irk<=Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos         
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    ds_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_ds

!!****f* input_cp2k_restarts/meta_hills_val_set_ww [1.0] *
!!
!!   NAME
!!     meta_hills_val_set_ww
!!
!!   FUNCTION
!!     routine to dump hills information during metadynamics run
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE meta_hills_val_set_ww(ww_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: ww_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_hills_val_set_ww', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(ww_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ww_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => ww_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(ww_section%values,2)==1) EXIT
       CALL section_vals_add_values(ww_section,error=error)
    END DO
    vals  => ww_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = meta_env%hills_env%n_hills
    DO irk=1,lsize
       CALL val_create(my_val,r_val=meta_env%hills_env%ww_history(irk),error=error)
       
       IF (irk<=Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos         
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    ww_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_ww

END MODULE input_cp2k_restarts



