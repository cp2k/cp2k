!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_electrostatic_methods [1.0] *
!!
!!   NAME MODULE
!!     qmmm_electrostatic_methods
!!
!!   FUNCTION
!!     A collection of methods to treat the QM/MM electrostatic coupling
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     5.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_electrostatic_methods
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dbl
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_init,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE qmmm_types,                      ONLY: qmmm_env_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE simulation_cell,                 ONLY: cell_type
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_electrostatic_methods'

  PUBLIC :: qmmm_el_coupling
!***
!****************************************************************************
CONTAINS
!!****f* qmmm_electrostatic_methods/qmmm_el_coupling [1.0] *
!!
!!   NAME
!!     qmmm_el_coupling
!!
!!   FUNCTION
!!     Main Driver to compute the QM/MM Electrostatic Coupling
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_el_coupling(qs_env,qmmm_env,particles_mm,error)
    Implicit None
    TYPE(qs_environment_type), INTENT(INOUT), &
         TARGET                                  :: qs_env   
    TYPE(qmmm_env_type), INTENT(in), &
         TARGET                                  :: qmmm_env  
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error
    TYPE(particle_type), DIMENSION(:), POINTER   :: particles_mm
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_el_coupling', &
         routineP = moduleN//':'//routineN    
    TYPE(cell_type), POINTER                 :: cell
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    LOGICAL :: failure
    INTEGER :: handle

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(cell, ks_qmmm_env_loc)
    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    error=error)
    ks_qmmm_env_loc => qs_env%ks_qmmm_env
    Write(6,'(T3,A)')"Information on the QM/MM Electrostatic Potential:"
    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE("GRID")
       Write(6,'(T3,A)')"QM/MM Coupling computed on GRID."
       ! Zeroing v_qmmm_rspace 
       CALL coeff_zero(coeff=ks_qmmm_env_loc%v_qmmm_rspace)
       ! Compute the MM electrostatic potential on GRID
       CALL qmmm_electrostatic_on_grid( v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                        mm_particles=particles_mm,&
                                        mm_charges=qmmm_env%mm_atom_chrg,&
                                        mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                        mm_atom_index=qmmm_env%mm_atom_index,&
                                        num_mm_atoms=qmmm_env%num_mm_atoms,&
                                        cell=cell,&
                                        Nexp=qmmm_env%grid_pot%Exponential_pot,&
                                        error=error)
       qs_env%ks_qmmm_env => ks_qmmm_env_loc
    CASE("GAUSS")
       Write(6,'(T3,A)')"QM/MM Coupling computed analytically within GTO."
       stop
    CASE("QS")
       Write(6,'(T3,A)')"QM/MM Coupling computed with the dual scheme GTO/GRID."
       stop
    CASE("NONE")
       Write(6,'(T3,A)')"No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
    CASE DEFAULT
       Write(6,'(T3,A)')"Unknown Coupling..."
       CPPrecondition(.false.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_el_coupling

!!****f* qmmm_electrostatic_methods/qmmm_electrostatic_on_grid [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_on_grid
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction on Grid Points
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_electrostatic_on_grid(v_qmmm,mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,cell,Nexp,error)
    USE  qmmm_coulomb_potential,  ONLY: Vpot
    Implicit None
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_electrostatic_on_grid', &
         routineP = moduleN//':'//routineN
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms, Nexp
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL :: failure
    INTEGER :: handle
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n
    INTEGER :: i, j, k, Imm, IndMM
    REAL(dbl) :: r, r0
    REAL(dbl), DIMENSION(3) :: xs, rv
   
    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)

    n(1:3) = v_qmmm%pw%pw_grid%npts_local (1:3)
    bo = v_qmmm%pw%pw_grid%bounds_local  
    
    MainLoopMM: Do Imm = 1, num_mm_atoms
       IndMM = mm_atom_index(Imm)
       r0 = mm_el_pot_radius(Imm) 
       !$omp parallel do private(k,j,i)
       LoopOnGrid: Do k = bo(1,3), bo(2,3)
          xs(3) = float(k)/float(n(3))
          Do j = bo(1,2), bo(2,2)
             xs(2) = float(j)/float(n(2))
             Do i = bo(1,1), bo(2,1)
                xs(1) =  float(i)/float(n(1))
                rv = mm_particles(IndMM)%r(:) - matmul(cell%hmat,xs) 
                r  = rv(1)*rv(1) + rv(2)*rv(2) + rv(3)*rv(3)
                v_qmmm%pw%cr3d(i,j,k) = v_qmmm%pw%cr3d(i,j,k) + &
                                        mm_charges(Imm) * Vpot(Nexp,r,r0)
             End Do
          End Do
       End Do LoopOnGrid
    End Do MainLoopMM
    

    CALL timestop(0.0_dp,handle)
  End SUBROUTINE qmmm_electrostatic_on_grid

End MODULE qmmm_electrostatic_methods
