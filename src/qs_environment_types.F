!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_environment_types [1.0] *
!!
!!   NAME
!!     qs_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (23.01.2002)
!!
!!   MODIFICATION HISTORY
!!     - mo_set_p_type added to qs_env (23.04.02,MK)
!!     - qs_force_type added to qs_env (05.06.02,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_environment_types

  USE kinds,                  ONLY: wp => dp

  USE atomic_kind_types,      ONLY: atomic_kind_type
  USE cp_b_matrix_structure,  ONLY: cp_b_matrix_struct_type,&
                                    cp_bmstruct_create,&
                                    cp_bmstruct_release
  USE cp_error_handling,      ONLY: cp_assert,&
                                    cp_assertion_failed,&
                                    cp_debug,&
                                    cp_error_dealloc_ref,&
                                    cp_error_get_logger,&
                                    cp_error_init,&
                                    cp_error_message,&
                                    cp_error_type,&
                                    cp_internal_error,&
                                    cp_unimplemented_error
  USE cp_fm_pool,             ONLY: cp_fm_pool_type,&
                                    fm_pool_create,&
                                    fm_pool_release,&
                                    fm_pool_retain
  USE cp_f_matrix_struct,     ONLY: cp_f_matrix_struct_type,&
                                    cp_fmstruct_create,&
                                    cp_fmstruct_release,&
                                    cp_fmstruct_retain
  USE cp_log_handling,        ONLY: cp_failure_level,&
                                    cp_fatal_level,&
                                    cp_log,&
                                    cp_logger_type,&
                                    cp_logger_get_default_unit_nr,&
                                    cp_note_level,&
                                    cp_to_string,&
                                    cp_warning_level
  USE cp_matrix_utils,        ONLY: cp_sm_set_dist_sparsity_from_sm
  USE cp_para_env,            ONLY: cp_para_env_release,&
                                    cp_para_env_retain
  USE cp_para_types,          ONLY: cp_blacs_env_type,&
                                    cp_para_env_type
  USE dft_types,              ONLY: dft_control_type
  USE pao_types,              ONLY: pao_env_type
  USE particle_types,         ONLY: particle_type
  USE pw_grid_types,          ONLY: HALFSPACE,&
                                    pw_grid_type
  USE pw_grids,               ONLY: pw_grid_construct,&
                                    pw_grid_setup
  USE pw_pool,                ONLY: pw_pool_create,&
                                    pw_pool_p_type,&
                                    pw_pool_release,&
                                    pw_pool_retain,&
                                    pw_pool_type,&
                                    pw_pools_copy,&
                                    pw_pools_dealloc
  USE qs_energy_types,        ONLY: qs_energy_type
  USE qs_force_types,         ONLY: qs_force_type
  USE qs_mo_types,            ONLY: mo_set_p_type
  USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
  USE qs_rho_types,           ONLY: qs_rho_type,&
                                    qs_rho_retain,&
                                    qs_rho_release
  USE scf_control_types,      ONLY: scf_control_type,&
                                    scf_c_release,&
                                    scf_c_retain
  USE simulation_cell,        ONLY: cell_type
  USE sparse_matrix_types,    ONLY: get_matrix_info,&
                                    real_matrix_p_type,&
                                    real_matrix_type
  USE timings,                ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_environment_types'

  TYPE qs_environment_type
!   PRIVATE
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(cell_type), POINTER                              :: cell
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: c
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_erfc,&
                                                             sac_ppl,sca_ppnl,&
                                                             scb_ppnl
    TYPE(pao_env_type), POINTER                           :: pao_env
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(qs_energy_type), POINTER                         :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,k,p,p_old,s,kinetic
    TYPE(cp_b_matrix_struct_type), POINTER                :: s_mstruct
    TYPE(qs_rho_type), POINTER                            :: rho
    TYPE(scf_control_type), POINTER                       :: scf_control

    TYPE(cp_fm_pool_type), POINTER                        :: mo_ao_fm_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER           :: pw_pools

    TYPE(cp_para_env_type), POINTER                       :: para_env
  END TYPE qs_environment_type

! *** Public data types ***

  PUBLIC :: qs_environment_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            init_qs_env,&
            set_qs_env,&
            qs_env_update_s_mstruct

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s, s_mstruct,&
                        rho, mo_ao_fm_pool, pw_pools, para_env, scf_control,&
                        kinetic,error)

!   Purpose: Get the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(IN)       :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER          :: cell
    TYPE(dft_control_type), OPTIONAL, POINTER   :: dft_control
    TYPE(mo_set_p_type),  DIMENSION(:),&
      OPTIONAL, POINTER                         :: c
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), OPTIONAL, POINTER           :: sab_orb,sac_erfc,sac_ppl,&
                                                   sca_ppnl,scb_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER       :: pao_env
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER     :: energy
    TYPE(qs_force_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: force
    TYPE(real_matrix_p_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: h,k,p,p_old,s,kinetic
    TYPE(cp_b_matrix_struct_type), OPTIONAL,&
         POINTER                                :: s_mstruct
    TYPE(qs_rho_type), OPTIONAL, POINTER        :: rho
    TYPE(cp_fm_pool_type), OPTIONAL, POINTER    :: mo_ao_fm_pool
    TYPE(pw_pool_p_type), OPTIONAL, &
         DIMENSION(:), POINTER                  :: pw_pools
    TYPE(cp_para_env_type), OPTIONAL, POINTER   :: para_env
    TYPE(scf_control_type), OPTIONAL, POINTER   :: scf_control
    TYPE(cp_error_type), OPTIONAL,INTENT(inout) :: error
!   ---------------------------------------------------------------------------

    IF (PRESENT(atomic_kind_set)) atomic_kind_set => qs_env%atomic_kind_set
    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(c)) c => qs_env%c
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sac_erfc)) sac_erfc => qs_env%sac_erfc
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sca_ppnl)) sca_ppnl => qs_env%sca_ppnl
    IF (PRESENT(scb_ppnl)) scb_ppnl => qs_env%scb_ppnl
    IF (PRESENT(pao_env)) pao_env => qs_env%pao_env
    IF (PRESENT(particle_set)) particle_set => qs_env%particle_set
    IF (PRESENT(energy)) energy => qs_env%energy
    IF (PRESENT(force)) force => qs_env%force
    IF (PRESENT(h)) h => qs_env%h
    IF (PRESENT(k)) k => qs_env%k
    IF (PRESENT(p)) p => qs_env%p
    IF (PRESENT(kinetic)) kinetic => qs_env%kinetic
    IF (PRESENT(p_old)) p_old => qs_env%p_old
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(s_mstruct)) s_mstruct => qs_env%s_mstruct
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(pw_pools)) pw_pools => qs_env%pw_pools
    IF (PRESENT(mo_ao_fm_pool)) mo_ao_fm_pool => qs_env%mo_ao_fm_pool
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env,para_env,error)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)       :: qs_env
    TYPE(cp_para_env_type), POINTER              :: para_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   ---------------------------------------------------------------------------

    NULLIFY (qs_env%atomic_kind_set)
    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%c)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sac_erfc)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sca_ppnl)
    NULLIFY (qs_env%scb_ppnl)
    NULLIFY (qs_env%pao_env)
    NULLIFY (qs_env%particle_set)
    NULLIFY (qs_env%energy)
    NULLIFY (qs_env%force)
    NULLIFY (qs_env%h)
    NULLIFY (qs_env%k)
    NULLIFY (qs_env%p)
    NULLIFY (qs_env%kinetic)
    NULLIFY (qs_env%p_old)
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%s_mstruct)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%mo_ao_fm_pool)
    NULLIFY (qs_env%pw_pools)
    NULLIFY (qs_env%scf_control)

    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env

  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s,&
                        pw_pools,mo_ao_fm_pool,rho,kinetic,&
                        scf_control,error)

!   Purpose: Set the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)    :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER        :: cell
    TYPE(dft_control_type), OPTIONAL, POINTER :: dft_control
    TYPE(mo_set_p_type),  DIMENSION(:),&
      OPTIONAL, POINTER                       :: c
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), OPTIONAL, POINTER         :: sab_orb,sac_erfc,sac_ppl,&
                                                 sca_ppnl,scb_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER     :: pao_env
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER   :: energy
    TYPE(qs_force_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: force
    TYPE(real_matrix_p_type), DIMENSION(:),&
         OPTIONAL, POINTER                    :: h,k,p,p_old,s,kinetic
    TYPE(qs_rho_type), POINTER,OPTIONAL       :: rho
    TYPE(cp_fm_pool_type), OPTIONAL, POINTER  :: mo_ao_fm_pool
    TYPE(pw_pool_p_type), OPTIONAL, &
         DIMENSION(:), POINTER                :: pw_pools
    TYPE(scf_control_type), OPTIONAL, POINTER :: scf_control
    TYPE(cp_error_type),INTENT(inout),OPTIONAL:: error

!   ---------------------------------------------------------------------------
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: new_pools
    NULLIFY(new_pools)

    IF (PRESENT(atomic_kind_set)) qs_env%atomic_kind_set => atomic_kind_set
    IF (PRESENT(cell)) qs_env%cell => cell
    IF (PRESENT(dft_control)) qs_env%dft_control => dft_control
    IF (PRESENT(c)) qs_env%c => c
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sac_erfc)) qs_env%sac_erfc => sac_erfc
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sca_ppnl)) qs_env%sca_ppnl => sca_ppnl
    IF (PRESENT(scb_ppnl)) qs_env%scb_ppnl => scb_ppnl
    IF (PRESENT(pao_env)) qs_env%pao_env => pao_env
    IF (PRESENT(particle_set)) qs_env%particle_set => particle_set
    IF (PRESENT(energy)) qs_env%energy => energy
    IF (PRESENT(force)) qs_env%force => force
    IF (PRESENT(h)) qs_env%h => h
    IF (PRESENT(k)) qs_env%k => k
    IF (PRESENT(p)) qs_env%p => p
    IF (PRESENT(p_old)) qs_env%p_old => p_old
    IF (PRESENT(s)) qs_env%s => s
    IF (PRESENT(rho)) THEN ! accept also null pointers?
      CALL qs_rho_retain(rho,error=error)
      CALL qs_rho_release(qs_env%rho,error=error)
      qs_env%rho => rho
    END IF
    IF (PRESENT(pw_pools)) THEN
      CALL pw_pools_copy(pw_pools,new_pools,error=error)
      CALL pw_pools_dealloc(qs_env%pw_pools,error=error)
      qs_env%pw_pools => new_pools
    END IF
    IF (PRESENT(mo_ao_fm_pool)) THEN ! accept also null pointers?
      CALL fm_pool_retain(mo_ao_fm_pool,error=error)
      CALL fm_pool_release(qs_env%mo_ao_fm_pool,error=error)
      qs_env%mo_ao_fm_pool => mo_ao_fm_pool
    END IF
    IF (PRESENT(scf_control)) THEN
       CALL scf_c_retain(scf_control,error=error)
       CALL scf_c_release(qs_env%scf_control,error=error)
       qs_env%scf_control => scf_control
    END IF
    IF (PRESENT(kinetic)) qs_env%kinetic => kinetic

  END SUBROUTINE set_qs_env

! *****************************************************************************

!!****f* qs_environment_types/qs_env_update_s_mstruct [1.0] *
!!
!!   NAME
!!     qs_env_update_s_mstruct
!!
!!   FUNCTION
!!     updates the s_mstruct to reflect the new overlap structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environment to update
!!     - para_env: the mpi environment
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_update_s_mstruct(qs_env,para_env,error)
  TYPE(qs_environment_type), INTENT(inout)     :: qs_env
  TYPE(cp_para_env_type), POINTER              :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_update_s_mstruct',&
        routineP=moduleN//':'//routineN
  INTEGER :: n_atoms, stat
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_sizes
  INTEGER, DIMENSION(:), POINTER :: first_row, first_col, last_row, last_col
  failure=.FALSE.

  IF (ASSOCIATED(qs_env%s_mstruct)) CALL cp_bmstruct_release(qs_env%s_mstruct)
  n_atoms=SIZE(qs_env%particle_set)
  ALLOCATE(block_sizes(n_atoms),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  block_sizes=qs_env%particle_set(:)%last_sgf-qs_env%particle_set(:)%first_sgf
  CALL cp_bmstruct_create(qs_env%s_mstruct,&
       symmetric=.TRUE.,&
       para_env=para_env,&
       b_row_sizes=block_sizes,&
       b_col_sizes=block_sizes,&
       error=error)
  CALL get_matrix_info(qs_env%s(1)%matrix,first_row=first_row,&
       last_row=last_row, first_col=first_col,&
       last_col=last_col)

  CALL cp_assert(ALL(qs_env%s_mstruct%b_dims%block_begins_at_row(1:n_atoms)==&
       first_row),cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  CALL cp_assert(ALL(qs_env%s_mstruct%b_dims%block_begins_at_col(1:n_atoms)==&
       first_col),cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  CALL cp_assert(qs_env%s_mstruct%b_dims%block_begins_at_col(n_atoms+1)==&
       last_col(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  CALL cp_assert(qs_env%s_mstruct%b_dims%block_begins_at_row(n_atoms+1)==&
       last_row(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)

  ! at the moment read from s, maybe in the future s will take its structure
  ! from this
  CALL cp_sm_set_dist_sparsity_from_sm(qs_env%s_mstruct,&
       para_env=para_env,&
       sample_matrix=qs_env%s(1)%matrix,&
       error=error)
END SUBROUTINE qs_env_update_s_mstruct
!***************************************************************************

!!****f* qs_environment_types/qs_env_rebuild_pw_pools [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_pw_pools
!!
!!   FUNCTION
!!     rebuilds the pw pools
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs environment where to build the pools
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_rebuild_pw_pools(qs_env, error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)     :: error

  LOGICAL :: failure, pw_grid_information
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_pw_pools',&
       routineP=moduleN//':'//routineN
  REAL(wp), DIMENSION(:), POINTER :: cutoff
  INTEGER :: handle,i,igrid_level,stat, ngrid_level
  TYPE(cell_type), POINTER          :: cell
  TYPE(pw_grid_type), POINTER :: pw_grid, old_pw_grid
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(pw_pool_p_type),DIMENSION(:), POINTER :: pw_pools
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_logger_type), pointer :: logger
  failure=.FALSE.
  NULLIFY(cutoff,cell,pw_grid,old_pw_grid,dft_control,pw_pools,para_env)

  CALL timeset(routineN//', '//moduleN,"I","",handle)
  CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control,&
       cell=cell)
  cutoff => dft_control%qs_control%e_cutoff
  ngrid_level=SIZE(cutoff)
  ALLOCATE(pw_pools(ngrid_level),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  logger=>cp_error_get_logger(error)
  pw_grid_information=logger%print_keys%pw_grid_information

  DO igrid_level=1,ngrid_level
     ALLOCATE(pw_grid,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     CALL pw_grid_construct(pw_grid)
     pw_grid%grid_span = HALFSPACE
     pw_grid%para%rs_dims ( 1 ) = para_env%num_pe
     pw_grid%para%rs_dims ( 2 ) = 1

     IF (igrid_level == 1) THEN
        IF (pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group,&
                info= cp_logger_get_default_unit_nr(&
                cp_error_get_logger(error),local=.FALSE.))
        ELSE
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group)
        END IF
     ELSE
        IF (pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group,&
                info=cp_logger_get_default_unit_nr(&
                cp_error_get_logger(error)),&
                ref_grid=old_pw_grid)
        ELSE
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group,&
                ref_grid=old_pw_grid)
        END IF
     END IF
     old_pw_grid => pw_grid

     NULLIFY(pw_pools(igrid_level)%pool)
     CALL pw_pool_create(pw_pools(igrid_level)%pool,&
          pw_grid_ptr=pw_grid,error=error)
  END DO

  CALL set_qs_env(qs_env, pw_pools=pw_pools)
  CALL pw_pools_dealloc(pw_pools,error=error)

  CALL timestop(0.0_wp,handle)
END SUBROUTINE qs_env_rebuild_pw_pools
!***************************************************************************

!!****f* qs_environment_types/qs_env_rebuild_fm_pools [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_fm_pools
!!
!!   FUNCTION
!!     rebuilds the pools of the full matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environement where the pools should be rebuilt
!!     - blacs_env: the blacs_env to use for the matrixes
!!     - para_env: a para_env coordinated to the blacs_env.
!!       Defaults to qs_env%para_env.
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_rebuild_fm_pools(qs_env, blacs_env, para_env,error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_blacs_env_type), POINTER                 :: blacs_env
  TYPE(cp_para_env_type), POINTER, OPTIONAL        :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)     :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_fm_pools',&
        routineP=moduleN//':'//routineN
  INTEGER :: handle
  TYPE(dft_control_type), POINTER         :: dft_control
  TYPE(cp_fm_pool_type), POINTER          :: moao_pool
  TYPE(cp_para_env_type), POINTER         :: p_env
  TYPE(cp_f_matrix_struct_type), POINTER  :: moao_mstruct
  failure=.FALSE.
  NULLIFY(dft_control, moao_pool)

  CALL timeset(routineN//', '//moduleN,"I","",handle)
  CALL get_qs_env(qs_env, para_env=p_env, dft_control=dft_control)
  IF (PRESENT(para_env)) p_env => para_env

!!FM  call cp_fmstruct_create(moao_mstruct, para_env=p_env,&
!!FM       context=blacs_env
!!FM       nrow_global=mo,ncol_global=ao)
!  call set_qs_env(qs_env, fm_pools=fm_pools)

  CALL timestop(0.0_wp,handle)
END SUBROUTINE qs_env_rebuild_fm_pools
!***************************************************************************

END MODULE qs_environment_types
