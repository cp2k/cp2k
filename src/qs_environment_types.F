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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_environment_types

  USE atomic_kind_types,      ONLY: atomic_kind_type
  USE dft_types,              ONLY: dft_control_type
  USE particle_types,         ONLY: particle_type
  USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
  USE simulation_cell,        ONLY: cell_type
  USE sparse_matrix_types,    ONLY: real_matrix_p_type

  IMPLICIT NONE

  PRIVATE

  TYPE qs_environment_type
    PRIVATE
    TYPE(cell_type), POINTER                              :: cell
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(real_matrix_p_type), POINTER                     :: h,s
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: k,p
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_erfc,&
                                                             sac_ppl,sca_ppnl,&
                                                             scb_ppnl
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

  SUBROUTINE get_qs_env(qs_env,cell,dft_control,h,k,p,s,atomic_kind_set,&
                        particle_set,sab_orb,sac_erfc,sac_ppl,sca_ppnl,&
                        scb_ppnl)

!   Purpose: Get the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(IN)         :: qs_env
    TYPE(cell_type), OPTIONAL, POINTER            :: cell
    TYPE(dft_control_type), OPTIONAL, POINTER     :: dft_control
    TYPE(real_matrix_p_type), OPTIONAL, POINTER   :: h,s
    TYPE(real_matrix_p_type), DIMENSION(:),&
      OPTIONAL, POINTER                           :: k,p
    TYPE(atomic_kind_type), DIMENSION(:),&
      OPTIONAL, POINTER                           :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                           :: particle_set
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), OPTIONAL, POINTER             :: sab_orb,sac_erfc,sac_ppl,&
                                                     sca_ppnl,scb_ppnl

!   ---------------------------------------------------------------------------

    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(h)) h => qs_env%h
    IF (PRESENT(k)) k => qs_env%k
    IF (PRESENT(p)) p => qs_env%p
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => qs_env%atomic_kind_set
    IF (PRESENT(particle_set)) particle_set => qs_env%particle_set
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sac_erfc)) sac_erfc => qs_env%sac_erfc
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sca_ppnl)) sca_ppnl => qs_env%sca_ppnl
    IF (PRESENT(scb_ppnl)) scb_ppnl => qs_env%scb_ppnl

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT) :: qs_env

!   ---------------------------------------------------------------------------

    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%h)
    NULLIFY (qs_env%k)
    NULLIFY (qs_env%p)
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%atomic_kind_set)
    NULLIFY (qs_env%particle_set)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sac_erfc)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sca_ppnl)
    NULLIFY (qs_env%scb_ppnl)

  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,cell,dft_control,h,k,p,s,atomic_kind_set,&
                        particle_set,sab_orb,sac_erfc,sac_ppl,sca_ppnl,&
                        scb_ppnl)

!   Purpose: Set the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)        :: qs_env
    TYPE(cell_type), OPTIONAL, POINTER            :: cell
    TYPE(dft_control_type), OPTIONAL, POINTER     :: dft_control
    TYPE(real_matrix_p_type), OPTIONAL, POINTER   :: h,s
    TYPE(real_matrix_p_type), DIMENSION(:),&
      OPTIONAL, POINTER                           :: k,p
    TYPE(atomic_kind_type), DIMENSION(:),&
      OPTIONAL, POINTER                           :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                           :: particle_set
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), OPTIONAL, POINTER             :: sab_orb,sac_erfc,sac_ppl,&
                                                     sca_ppnl,scb_ppnl

!   ---------------------------------------------------------------------------

    IF (PRESENT(cell)) qs_env%cell => cell
    IF (PRESENT(dft_control)) qs_env%dft_control => dft_control
    IF (PRESENT(h)) qs_env%h => h
    IF (PRESENT(k)) qs_env%k => k
    IF (PRESENT(p)) qs_env%p => p
    IF (PRESENT(s)) qs_env%s => s
    IF (PRESENT(atomic_kind_set)) qs_env%atomic_kind_set => atomic_kind_set
    IF (PRESENT(particle_set)) qs_env%particle_set => particle_set
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sac_erfc)) qs_env%sac_erfc => sac_erfc
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sca_ppnl)) qs_env%sca_ppnl => sca_ppnl
    IF (PRESENT(scb_ppnl)) qs_env%scb_ppnl => scb_ppnl

  END SUBROUTINE set_qs_env

! *****************************************************************************

END MODULE qs_environment_types
