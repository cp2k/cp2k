MODULE scp_restarts

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: energy_force_run,&
                                             energy_run
  USE input_cp2k_restarts,             ONLY: write_restart_header
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_restarts'
  
  PUBLIC :: update_scp_section, write_scp_restart

  CONTAINS

! *****************************************************************************
  SUBROUTINE write_scp_restart (scp_env, root_section, globenv, atomic_kind_set, particle_set, error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_scp_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=30), DIMENSION(2) :: &
      keys = (/"SCF%PRINT%RESTART_HISTORY","SCF%PRINT%RESTART        "/)
    INTEGER                                  :: ikey, ires
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: dft_section

       logger => cp_error_get_logger(error)
       CALL get_scp_env(scp_env=scp_env,&
            input=dft_section,&
            error=error)

! Test against the DFT%SCF print keys for writing wavefunctions.
       IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(1),error=error),cp_p_file) .OR.  &
         BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(2),error=error),cp_p_file) ) THEN
            
         CALL update_scp_section ( scp_env, atomic_kind_set, particle_set, error ) 
         IF (globenv%run_type_id == energy_run .OR. globenv%run_type_id==energy_force_run) THEN
            DO ikey=1,SIZE(keys)
              IF (BTEST(cp_print_key_should_output(logger%iter_info,&
                 dft_section,keys(ikey),error=error),cp_p_file)) THEN

                 ires = cp_print_key_unit_nr(logger,dft_section,TRIM(keys(ikey)),&
                        extension=".restart",  do_backup=(ikey==2), error=error)
                 IF (ires>0) THEN
                    CALL write_restart_header(ires, error)
                    CALL section_vals_write(root_section,unit_nr=ires,hide_root=.FALSE., error=error)
                 ENDIF

                 CALL cp_print_key_finished_output(ires,logger,dft_section,&
                      TRIM(keys(ikey)), error=error)
              END IF
            END DO
         END IF
       END IF
  END SUBROUTINE write_scp_restart 
      
! *****************************************************************************
  SUBROUTINE update_scp_section (scp_env, atomic_kind_set, particle_set,  error)

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_scp_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, icoef, iel, &
                                                ikind, irep, n_el, natom, &
                                                natom_tot, nkind, nsgf, &
                                                nsgf_max, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cindex
    INTEGER, POINTER                         :: INDEX( : )
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbuffer, pindex
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: coeffs
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(section_vals_type), POINTER         :: dft_section, scp_section, &
                                                work_section

       CALL timeset(routineN,"I","",handle)

       failure = .FALSE.
       logger => cp_error_get_logger(error)
       NULLIFY ( dft_section, scp_section, para_env, INDEX )

       CALL get_scp_env(scp_env=scp_env,&
            aux_coeff_set=aux_coeff_set,&
            input=dft_section,&
            para_env=para_env,&
            error=error)

       
       IF (ASSOCIATED(aux_coeff_set)) THEN          
          natom_tot=SIZE ( particle_set )
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
          ALLOCATE (cindex(natom_tot),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (INDEX(n_el),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (pindex(natom_tot,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   

          pindex(:,:)=0.0_dp
          cindex(:)=0

          nkind=SIZE(aux_coeff_set%coeffs_of_kind)
          DO ikind=1, nkind
             coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
             atomic_kind=> atomic_kind_set(ikind)
             IF(.NOT.ASSOCIATED(coeffs)) CYCLE
             DO iel= 1,coeffs%n_els
                iatom=aux_coeff_set%distribution%list(ikind)%array(iel)
                cindex(iatom)=iatom
                DO icoef=1,coeffs%ncoef_atom
                    pindex(iatom,icoef)=coeffs%c(iel,icoef)
                END DO
             END DO
          END DO
          CALL mp_sum(pindex,para_env%group)
          CALL mp_sum(cindex, para_env%group)
        
          iel = 0
          DO iatom = 1, natom_tot
            IF ( cindex ( iatom ) == 0 ) THEN
              CYCLE 
            ELSE
              iel = iel + 1
              INDEX(iel) = cindex ( iatom )
              pbuffer(iel,:) = pindex ( iatom, : )
            END IF
          END DO

          scp_section => section_vals_get_subs_vals(dft_section,"SCP",error=error)
          work_section => section_vals_get_subs_vals(scp_section,"SCP_COEFF%RESTART_INDEX",error=error)
          CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",i_vals_ptr=index,error=error)
          work_section => section_vals_get_subs_vals(scp_section,"SCP_COEFF%COORD",error=error)
          DO irep = 1, SIZE(index)
             ALLOCATE(buffer(SIZE(pbuffer,2)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = pbuffer(irep,:)
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_vals_ptr=buffer,&
                  i_rep_val=irep,error=error)    
          END DO
          DEALLOCATE (pindex,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          DEALLOCATE (cindex,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (pbuffer,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
       END IF
  
       CALL timestop(0.0_dp,handle)

  END SUBROUTINE update_scp_section

! *****************************************************************************
END MODULE scp_restarts
