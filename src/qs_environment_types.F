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

  USE atomic_kind_types,      ONLY: atomic_kind_type
  USE dft_types,              ONLY: dft_control_type
  USE pao_types,              ONLY: pao_env_type
  USE particle_types,         ONLY: particle_type
  USE qs_energy_types,        ONLY: qs_energy_type
  USE qs_force_types,         ONLY: qs_force_type
  USE qs_mo_types,            ONLY: mo_set_p_type
  USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
  USE simulation_cell,        ONLY: cell_type
  USE sparse_matrix_types,    ONLY: real_matrix_p_type, real_matrix_type, get_matrix_info
  USE cp_b_matrix_structure,  ONLY: cp_b_matrix_struct_type,&
       cp_bmstruct_create, cp_bmstruct_release
  use cp_log_handling,        only: cp_to_string, cp_failure_level,&
       cp_warning_level,cp_note_level
  use cp_error_handling,      only: cp_assert, cp_error_type, &
       cp_assertion_failed
  use cp_matrix_utils,        only: cp_sm_set_dist_sparsity_from_sm
  use cp_para_env,            only: cp_para_env_type
  use qs_rho_types,           only: qs_rho_type
       

  IMPLICIT NONE

  PRIVATE

  character(len=*), private, parameter :: moduleN='qs_environment_types'

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
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,k,p,p_old,s
    TYPE(cp_b_matrix_struct_type), POINTER                :: s_mstruct
    TYPE(qs_rho_type), pointer                            :: rho
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
                        rho)

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
      OPTIONAL, POINTER                         :: h,k,p,p_old,s
    TYPE(cp_b_matrix_struct_type), OPTIONAL,&
         POINTER                                :: s_mstruct
    TYPE(qs_rho_type), optional, pointer       :: rho

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
    IF (PRESENT(p_old)) p_old => qs_env%p_old
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(s_mstruct)) s_mstruct => qs_env%s_mstruct
    IF (PRESENT(rho)) rho => qs_env%rho

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT) :: qs_env

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
    NULLIFY (qs_env%p_old)
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%s_mstruct)
    NULLIFY (qs_env%rho)

  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s)

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
      OPTIONAL, POINTER                       :: h,k,p,p_old,s

!   ---------------------------------------------------------------------------

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
subroutine qs_env_update_s_mstruct(qs_env,para_env,error)
  type(qs_environment_type), intent(inout) :: qs_env
  type(cp_para_env_type), pointer :: para_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_env_update_s_mstruct',&
        routineP=moduleN//':'//routineN
  integer :: n_atoms, stat
  integer, dimension(:), allocatable :: block_sizes
  integer, dimension(:), pointer :: first_row, first_col, last_row, last_col
  failure=.false.
  
  if (associated(qs_env%s_mstruct)) call cp_bmstruct_release(qs_env%s_mstruct)
  n_atoms=size(qs_env%particle_set)
  allocate(block_sizes(n_atoms),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  block_sizes=qs_env%particle_set(:)%last_sgf-qs_env%particle_set(:)%first_sgf
  call cp_bmstruct_create(qs_env%s_mstruct,&
       symmetric=.true.,&
       para_env=para_env,&
       b_row_sizes=block_sizes,&
       b_col_sizes=block_sizes,&
       error=error)
  call get_matrix_info(qs_env%s(1)%matrix,first_row=first_row,&
       last_row=last_row, first_col=first_col,&
       last_col=last_col)

  call cp_assert(all(qs_env%s_mstruct%b_dims%block_begins_at_row(1:n_atoms)==&
       first_row),cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  call cp_assert(all(qs_env%s_mstruct%b_dims%block_begins_at_col(1:n_atoms)==&
       first_col),cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  call cp_assert(qs_env%s_mstruct%b_dims%block_begins_at_col(n_atoms+1)==&
       last_col(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  call cp_assert(qs_env%s_mstruct%b_dims%block_begins_at_row(n_atoms+1)==&
       last_row(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)

  ! at the moment read from s, maybe in the future s will take its structure
  ! from this
  call cp_sm_set_dist_sparsity_from_sm(qs_env%s_mstruct,&
       para_env=para_env,&
       sample_matrix=qs_env%s(1)%matrix,&
       error=error)
end subroutine qs_env_update_s_mstruct
!***************************************************************************

END MODULE qs_environment_types
