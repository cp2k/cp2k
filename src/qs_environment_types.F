!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
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
  USE sparse_matrix_types,    ONLY: real_matrix_p_type

  IMPLICIT NONE

  PRIVATE

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
!MK should be removed:
    INTEGER                                               :: n_electrons
  END TYPE qs_environment_type

! *** Public data types ***

  PUBLIC :: qs_environment_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            init_qs_env,&
            set_qs_env

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s,&
                        natoms,n_electrons)

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
    INTEGER, INTENT(OUT), OPTIONAL              :: natoms,n_electrons

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
!MK should be retrieved directly
    IF (PRESENT(natoms)) natoms = SIZE(qs_env%particle_set)
    IF (PRESENT(n_electrons)) n_electrons = qs_env%n_electrons

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

  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s,&
                        n_electrons)

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
    INTEGER, INTENT(IN), OPTIONAL             :: n_electrons

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
    IF (PRESENT(n_electrons)) qs_env%n_electrons = n_electrons

  END SUBROUTINE set_qs_env

! *****************************************************************************

END MODULE qs_environment_types
