
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_methods [1.0] *
!!
!!   NAME
!!     qmmm_methods
!!
!!   FUNCTION
!!     Routines to compute energy and forces in a QM/MM calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!
!!   SOURCE
!****************************************************************************

MODULE qmmm_methods
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
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE memory_utilities,                ONLY: reallocate
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_init,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE qmmm_types,                      ONLY: qmmm_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_create_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_type,&
                                             pw_sumup,&
                                             pw_transfer
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE pws,                             ONLY: pw_integral_ab
  USE mathconstants,                   ONLY: pi
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_get_loop_vars,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_transfer,&
                                             rs_inject_NoPBC
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type
  USE util,                            ONLY: get_limit
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE mm_collocate_potential,          ONLY: integrate_gf_rspace_NoPBC

  IMPLICIT NONE
  
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_methods'
  PUBLIC :: qmmm_calculate_energy,&
            qmmm_modify_hartree_pot,&
            qmmm_forces
  
CONTAINS
!!****f* qmmm_methods/qmmm_calculate_energy [1.0] *
!!
!!   NAME
!!     qmmm_calculate_energy
!!
!!   FUNCTION
!!     Computes the contribution to the total energy of the QM/MM 
!!     electrostatic coupling
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_calculate_energy(qs_env, rho, v_qmmm, qmmm_energy, error)
    IMPLICIT NONE
    TYPE(qs_environment_type), INTENT(in), &
         TARGET                                 :: qs_env
    TYPE(coeff_type), DIMENSION(:), POINTER :: rho
    TYPE(coeff_type), INTENT(IN) :: v_qmmm
    REAL(dbl), INTENT(INOUT) :: qmmm_energy
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_calculate_energy', &
         routineP = moduleN//':'//routineN   
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL :: failure
    INTEGER :: handle, ispin
   
    CALL timeset(routineN,"I","",handle)    
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    NULLIFY(dft_control)

    call get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    error=error)

    qmmm_energy = 0.d0
    Do ispin=1, dft_control%nspins
       qmmm_energy = qmmm_energy + pw_integral_ab ( rho(ispin)%pw, v_qmmm%pw )
    End Do

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_calculate_energy

!!****f* qmmm_methods/qmmm_modify_hartree_pot [1.0] *
!!
!!   NAME
!!     qmmm_modify_hartree_pot
!!
!!   FUNCTION
!!     Modify the hartree potential in order to include the QM/MM correction
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_modify_hartree_pot(v_hartree,v_qmmm, error)
    IMPLICIT NONE
    TYPE(coeff_type), INTENT(INOUT) :: v_hartree
    TYPE(coeff_type), INTENT(IN)    :: v_qmmm
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_modify_hartree_pot', &
         routineP = moduleN//':'//routineN   
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    LOGICAL :: failure
    INTEGER :: handle
   
    CALL timeset(routineN,"I","",handle)    
    failure=.FALSE.

    v_hartree%pw%cr3d = v_hartree%pw%cr3d + &
                        v_qmmm%pw%pw_grid%dvol*v_qmmm%pw%cr3d

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_modify_hartree_pot

!!****f* qmmm_methods/qmmm_forces [1.0] *
!!
!!   NAME
!!     qmmm_forces
!!
!!   FUNCTION
!!     General driver to Compute the contribution 
!!     to the forces due to the QM/MM potential
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces(qs_env,qmmm_env,mm_particles,error)
    IMPLICIT NONE
    TYPE(qs_environment_type), INTENT(INOUT), &
         TARGET                                  :: qs_env   
    TYPE(qmmm_env_type), INTENT(in), &
         TARGET                                  :: qmmm_env  
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                   :: error
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
   
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces', &
         routineP = moduleN//':'//routineN    
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(pw_env_type),POINTER                :: pw_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
         POINTER                      :: rs_pools
    TYPE(pw_type),POINTER                    :: rho_tot_r
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_pool_type),POINTER               :: auxbas_pool
    LOGICAL :: failure
    INTEGER :: handle, ispin, auxbas_grid, Imm, IndMM

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(ks_qmmm_env_loc, rho, pw_env, rs_pools, rho_tot_r, energy)
    CALL get_qs_env(qs_env=qs_env,&
                    rho=rho,&
                    rho_core=rho_core,&
                    pw_env=pw_env,&
                    energy=energy,&
                    para_env=para_env,&
                    error=error)
    ks_qmmm_env_loc => qs_env%ks_qmmm_env
    !
    CALL pw_env_get(pw_env=pw_env,&
                    rs_pools=rs_pools,& 
                    auxbas_grid=auxbas_grid,&
                    auxbas_pw_pool=auxbas_pool,&
                    error=error)
    CALL pw_pool_create_pw(auxbas_pool,rho_tot_r,&
         in_space=REALSPACE,use_data=REALDATA3D, error=error)
    CALL pw_transfer(rho_core%pw,rho_tot_r)
    !
    ! Computes the QM/MM Nuclear Electrostatic Potential
    !
    energy%qmmm_nu = pw_integral_ab ( rho_tot_r, ks_qmmm_env_loc%v_qmmm_rspace%pw )
    energy%total = energy%total +  energy%qmmm_nu
    WRITE(6,'(T3,A,F15.9)')" QM/MM Nuclear Electrostatic Potential :",energy%qmmm_nu
    WRITE(6,'(T3,A,F15.9)')" Total Energy after the QM/MM Calculation :",energy%total    
    !
    DO ispin=1,SIZE(rho%rho_r)
       CALL pw_sumup(rho%rho_r(ispin)%pw,rho_tot_r)
    END DO
    WRITE(6,'(T3,A)')" Evaluating forces on MM atoms due to the:"
    !
    ! Electrostatic Interaction type...
    !
    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE("GRID")
       WRITE(6,'(T3,A)')"QM/MM Coupling computed on GRID."
       ! Compute the MM electrostatic potential on GRID
       CALL qmmm_forces_on_grid( rho=rho_tot_r,&
                                 v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                 mm_particles=mm_particles,&
                                 mm_charges=qmmm_env%mm_atom_chrg,&
                                 mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                 mm_atom_index=qmmm_env%mm_atom_index,&
                                 num_mm_atoms=qmmm_env%num_mm_atoms,&
                                 dvol=rho%rho_r(1)%pw%pw_grid%dvol,&
                                 error=error)
    CASE("GAUSS")
       Write(6,'(T3,A)')"QM/MM Coupling computed collocating the Gaussian Potential Functions."
       CALL qmmm_forces_with_gaussian(rho=rho_tot_r,&
                                      mm_particles=mm_particles,&
                                      mm_charges=qmmm_env%mm_atom_chrg,&
                                      mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                      mm_atom_index=qmmm_env%mm_atom_index,&
                                      num_mm_atoms=qmmm_env%num_mm_atoms,&
                                      rs_pools=rs_pools,&
                                      auxbas_grid=auxbas_grid,&
                                      para_env=para_env,&
                                      eps_mm_rspace=qmmm_env%eps_mm_rspace,&
                                      cube_info=ks_qmmm_env_loc%cube_info,&
                                      pgfs=qmmm_env%pgfs,&
                                      error=error)
    CASE("NONE")
       Write(6,'(T3,A)')"No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
    CASE DEFAULT
       Write(6,'(T3,A)')"Unknown Coupling..."
       CPPrecondition(.false.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_r,&
         accept_non_compatible=.TRUE.,error=error)

    ! Print Out Forces
    DO Imm = 1, SIZE(qmmm_env%mm_atom_index)
       IndMM = qmmm_env%mm_atom_index(Imm)
       WRITE(6,'(I6,3F15.9)') IndMM, mm_particles(IndMM)%f
    END DO
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_forces

!!****f* qmmm_methods/qmmm_forces_on_grid [1.0] *
!!
!!   NAME
!!     qmmm_forces_on_grid
!!
!!   FUNCTION
!!     Evaluates the contribution 
!!     to the forces due to the QM/MM potential computed on GRID
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_on_grid(rho,v_qmmm,mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index, num_mm_atoms, dvol,error)
    Implicit None
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_on_grid', &
         routineP = moduleN//':'//routineN
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    TYPE(pw_type),POINTER                        :: rho
    REAL(dbl), INTENT(IN) :: dvol
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms
    INTEGER :: handle
    
    CALL timeset(routineN,"I","",handle)
    CALL qmmm_forces_on_grid_low(v_qmmm,&
                                 rho%cr3d,&
                                 mm_particles,&
                                 mm_charges,&
                                 mm_el_pot_radius,&
                                 mm_atom_index,&
                                 num_mm_atoms,&
                                 dvol,&
                                 error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_forces_on_grid

!!****f* qmmm_methods/qmmm_forces_on_grid_low [1.0] *
!!
!!   NAME
!!     qmmm_forces_on_grid_low
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the 
!!     QM/MM potential computed on GRID. Low Level Routine.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_on_grid_low(v_qmmm,Rho,mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,dvol,error)
    Implicit None
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_on_grid_low', &
         routineP = moduleN//':'//routineN
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    REAL(dbl), DIMENSION(*) :: Rho
    REAL(dbl), INTENT(IN) :: dvol
    INTEGER, INTENT(IN) :: num_mm_atoms
    LOGICAL :: failure
    INTEGER :: handle
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER                   :: n1,n2,n3
    INTEGER :: i, j, k, Imm, IndMM,jj,my_k,my_j
    REAL(dbl) :: r, r0, qt, rNexp, rNexp0, rNexp1, r0Nexp, r0Nexp1, dpot
    REAL(dbl) :: dr1,dr2,dr3, xs1,xs2,xs3, rv1,rv2,rv3, rt1,rt2,rt3
    REAL(dbl) :: ft1, ft2, ft3, fac

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)

    bo = v_qmmm%pw%pw_grid%bounds_local  
    n1 = v_qmmm%pw%pw_grid%npts(1)   
    n2 = v_qmmm%pw%pw_grid%npts(2)   
    n3 = v_qmmm%pw%pw_grid%npts(3)   
    dr1 = v_qmmm%pw%pw_grid%dr(1)
    dr2 = v_qmmm%pw%pw_grid%dr(2)
    dr3 = v_qmmm%pw%pw_grid%dr(3)
    CALL timeset(routineN,"I","",handle)

    MainLoopMM: Do Imm = 1, num_mm_atoms
       IndMM = mm_atom_index(Imm)
       r0 = mm_el_pot_radius(Imm) 
       r0Nexp = r0**2
       r0Nexp = r0Nexp*r0Nexp
       r0Nexp1 = r0*r0Nexp
       rt1= mm_particles(IndMM)%r(1) 
       rt2= mm_particles(IndMM)%r(2) 
       rt3= mm_particles(IndMM)%r(3)
       ft1 =0.0_dbl
       ft2 =0.0_dbl
       ft3 =0.0_dbl
       qt=mm_charges(Imm)*dvol
       !$omp parallel do default (none) private( i, j, k, Imm, IndMM,jj,my_k,my_j,&
       !$omp                            r, r0, qt, rNexp, rNexp0, rNexp1, r0Nexp, r0Nexp1,&
       !$omp                   dr1,dr2,dr3, xs1,xs2,xs3, rv1,rv2,rv3, rt1,rt2,rt3,&
       !$omp                   ft1,ft2,ft3,dpot,fac) shared(Rho,bo,n1,n2,n3)
       LoopOnGrid: Do k = bo(1,3), bo(2,3)
          my_k=(k+n3/2)
          xs3 = real(my_k,dbl)*dr3
          my_j=bo(1,2)+n2/2
          jj=my_k*n1*n2+my_j*n1
          xs2 = REAL(my_j,dbl)*dr2
          rv3 = rt3 - xs3
          Do j = bo(1,2), bo(2,2)
             xs1 = (bo(1,1)+n1/2)*dr1
             rv2 = rt2 - xs2
             Do i = bo(1,1), bo(2,1)
                jj=jj+1
                rv1 = rt1 - xs1
                r  = sqrt(rv1*rv1 + rv2*rv2 + rv3*rv3)
                rNexp  = r**2
                rNexp0 = rNexp*r
                rNexp  = rNexp*rNexp
                rNexp1 = rNexp*r
                dpot = 4.0_dbl * (rNexp1-r0Nexp1)* rNexp0 -&
                       5.0_dbl * (rNexp -r0Nexp) * rNexp
                dpot = dpot / ((rNexp1 - r0Nexp1)**2 * r) 
                fac= - Rho(jj)* dpot 
                ft1 = ft1 + fac * rv1
                ft2 = ft2 + fac * rv2
                ft3 = ft3 + fac * rv3
                xs1 = xs1 + dr1
             End Do
             xs2 = xs2 + dr2
          End Do
       END DO LoopOnGrid
       mm_particles(IndMM)%f(1) = ft1 * qt + mm_particles(IndMM)%f(1) 
       mm_particles(IndMM)%f(2) = ft2 * qt + mm_particles(IndMM)%f(2) 
       mm_particles(IndMM)%f(3) = ft3 * qt + mm_particles(IndMM)%f(3) 
    END DO MainLoopMM

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_forces_on_grid_low
!!****f* qmmm_methods/qmmm_forces_with_gaussian [1.0] *
!!
!!   NAME
!!     qmmm_forces_with_gaussian
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the 
!!     QM/MM potential computed collocating the Electrostatic
!!     Gaussian Potential.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian(rho, mm_particles,&
       mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
       rs_pools, auxbas_grid, cube_info, para_env, eps_mm_rspace, pgfs,&
       error)
    IMPLICIT NONE
     ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(pw_type),POINTER                        :: rho
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms, auxbas_grid
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
         POINTER                                 :: rs_pools  
    REAL(KIND=dbl), INTENT(IN) :: eps_mm_rspace
    TYPE(cp_para_env_type), POINTER              :: para_env
    TYPE(cube_info_type), POINTER, DIMENSION(:)  :: cube_info
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                   :: pgfs
    ! Local Variables
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian', &
         routineP = moduleN//':'//routineN
    TYPE(realspace_grid_p_type), POINTER, DIMENSION(:) :: tmp_grid
    REAL(KIND=dp) :: dvol
    INTEGER :: ngrids, igrid, stat, handle
    ! Statements
    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    NULLIFY(tmp_grid)
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)
    !Statements
    ngrids=size(rs_pools)
    ALLOCATE(tmp_grid(ngrids),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO igrid=1,ngrids
       NULLIFY(tmp_grid(igrid)%rs_grid)
       CALL rs_pool_create_rs(rs_pools(igrid)%pool,tmp_grid(igrid)%rs_grid)
       CALL rs_grid_zero(tmp_grid(igrid)%rs_grid)
    END DO

    ! Collocate Density on multigrids
    CALL rs_pw_transfer(tmp_grid(auxbas_grid)%rs_grid,rho,"BACKWARD")
    DO Igrid = auxbas_grid, SIZE(tmp_grid)-1
       CALL rs_inject_NoPBC(tmp_grid(Igrid  )%rs_grid,&
                            tmp_grid(igrid+1)%rs_grid,&
                            error)
    END DO
    dvol = rho%pw_grid%dvol

    CALL qmmm_forces_with_gaussian_low( tmp_grid, mm_particles,&
         mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
         cube_info, para_env, eps_mm_rspace, dvol, pgfs, error)

    DO igrid=1,ngrids
       CALL rs_pool_give_back_rs(rs_pools(igrid)%pool,tmp_grid(igrid)%rs_grid,&
            error=error)
    END DO
    DEALLOCATE(tmp_grid,stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_forces_with_gaussian

!!****f* qmmm_methods/qmmm_forces_with_gaussian_low [1.0] *
!!
!!   NAME
!!     qmmm_forces_with_gaussian_low
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the 
!!     QM/MM potential computed collocating the Electrostatic
!!     Gaussian Potential. Low Level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian_low(tmp_grid, mm_particles,&
       mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
       cube_info, para_env, eps_mm_rspace, dvol, pgfs, error)
    IMPLICIT NONE
     ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(realspace_grid_p_type),POINTER, DIMENSION(:) :: tmp_grid
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms
    REAL(KIND=dbl), INTENT(IN) :: eps_mm_rspace, dvol
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), POINTER, DIMENSION(:)  :: cube_info
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                   :: pgfs
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian_low', &
         routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(3)              :: force
    LOGICAL :: failure
    INTEGER :: stat, handle, ithread, igrid, ilevel
    INTEGER, DIMENSION(3) :: ndim
    INTEGER :: Imm, IndMM, LIndMM, IRadTyp, Icg, IGauss
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp) :: alpha, height, W
    TYPE(qmmm_gaussian_type), POINTER :: pgf
    REAL(KIND=dbl), DIMENSION(:), POINTER    :: xdat, ydat, zdat
    REAL(KIND=dbl), DIMENSION(:), POINTER    :: xder, yder, zder, wrk
    !Statements
    NULLIFY(pgf, xdat, ydat, zdat, xder, yder, zder, wrk)    
    CALL timeset(routineN,"I","",handle)
    ithread = 0
    Icg     = 0
    ndim    = 0
    DO igrid=1,SIZE(tmp_grid)
       ndim = MAX(tmp_grid(igrid)%rs_grid%npts,ndim)
    END DO
    ALLOCATE (xdat(ndim(1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(ndim(2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(ndim(3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (xder(ndim(1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (yder(ndim(2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zder(ndim(3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (wrk(MAX(ndim(1),ndim(2),ndim(3))), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!
    DO IRadTyp = 1, SIZE(pgfs) 
       pgf => pgfs(IRadTyp)%pgf
!       WRITE(*,'(A,F12.6)')' Collocating MM Potential with Radius: ', pgf%My_Elp_Radius
       DO IGauss = 1, pgf%Number_of_Gaussians
          alpha     = 1.0_dbl / pgf%Gk(IGauss)
          alpha     = alpha * alpha
          height    = pgf%Ak(IGauss)
!          WRITE(*,'(5X,A,I5)')' Collocating Gaussian Number:',IGauss
          DO Imm = 1, SIZE(pgf%mm_atom_index)
             Icg       =   Icg + 1
             LIndMM    =   pgf%mm_atom_index(Imm)
             IndMM     =   mm_atom_index(LIndMM)
             ra(:)     =   mm_particles(IndMM)%r
             W         =   mm_charges(LIndMM) * height
             ilevel    =   pgf%grid_level(IGauss)
             force     =   0.0_dbl
!             WRITE(*,'(10X,I5,A,I5,A,3F12.6,A,F12.6,A,F12.6)')&
!                  Icg,' MM Atom:',IndMM,&
!                  ' x,y,z:',ra(1), ra(2), ra(3),&
!                  ' alpha:',alpha,&
!                  ' prefactor:',W

             CALL integrate_gf_rspace_NoPBC(zetp=alpha,&
                                            rp=ra,&
                                            scale=1.0_dp,&
                                            W=W,&
                                            rsgrid=tmp_grid(ilevel)%rs_grid,&
                                            cube_info=cube_info(ilevel),&
                                            eps_mm_rspace=eps_mm_rspace,&
                                            ithread=ithread, &
                                            xdat=xdat,&
                                            ydat=ydat,&
                                            zdat=zdat,&
                                            xder=xder,&
                                            yder=yder,&
                                            zder=zder,&
                                            wrk=wrk,&
                                            force=force)

             mm_particles(IndMM)%f(:) = force(:) * dvol + mm_particles(IndMM)%f(:) 
          END DO
       END DO
    END DO

    IF (ASSOCIATED(xdat)) THEN
       DEALLOCATE (xdat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(ydat)) THEN
       DEALLOCATE (ydat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(zdat)) THEN
       DEALLOCATE (zdat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(xder)) THEN
       DEALLOCATE (xder, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(yder)) THEN
       DEALLOCATE (yder, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(zder)) THEN
       DEALLOCATE (zder, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(wrk)) THEN
       DEALLOCATE (wrk, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_forces_with_gaussian_low

END MODULE qmmm_methods
