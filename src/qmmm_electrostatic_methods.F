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
  USE memory_utilities,                ONLY: reallocate
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
  USE mathconstants,                   ONLY: pi
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_get_loop_vars,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type
  USE l_utils,                         ONLY: l_info_type
  USE util,                            ONLY: get_limit
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace

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
    ! Arguments
    TYPE(qs_environment_type), INTENT(INOUT), &
         TARGET                                  :: qs_env   
    TYPE(qmmm_env_type), INTENT(in), &
         TARGET                                  :: qmmm_env  
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error
    TYPE(particle_type), DIMENSION(:), POINTER   :: particles_mm
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_el_coupling', &
         routineP = moduleN//':'//routineN    
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    LOGICAL :: failure
    INTEGER :: handle
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    ! Statements
    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(ks_qmmm_env_loc)
    CALL get_qs_env(qs_env=qs_env,&
                    pw_env=pw_env,&
                    para_env=para_env,&
                    error=error)
    CALL pw_env_get(pw_env,auxbas_rs_pool=auxbas_rs_pool,&
         error=error)
    ks_qmmm_env_loc => qs_env%ks_qmmm_env
    Write(6,'(T3,A)')"Information on the QM/MM Electrostatic Potential:"
    !
    ! Initializing vectors:
    !        Zeroing v_qmmm_rspace 
    CALL coeff_zero(coeff=ks_qmmm_env_loc%v_qmmm_rspace)
    !
    !
    ! Electrostatic Interaction type...
    !
    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE("GRID")
       Write(6,'(T3,A)')"QM/MM Coupling computed on GRID."
       ! Compute the MM electrostatic potential on GRID
       CALL qmmm_electrostatic_on_grid( v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                        mm_particles=particles_mm,&
                                        mm_charges=qmmm_env%mm_atom_chrg,&
                                        mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                        mm_atom_index=qmmm_env%mm_atom_index,&
                                        num_mm_atoms=qmmm_env%num_mm_atoms,&
                                        error=error)
       qs_env%ks_qmmm_env => ks_qmmm_env_loc
    CASE("GAUSS")
       Write(6,'(T3,A)')"QM/MM Coupling computed collocating the Gaussian Potential Functions."
       CALL    qmmm_elec_with_gaussian(v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                              mm_particles=particles_mm,&
                                              mm_charges=qmmm_env%mm_atom_chrg,&
                                              mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                              mm_atom_index=qmmm_env%mm_atom_index,&
                                              num_mm_atoms=qmmm_env%num_mm_atoms,&
                                              auxbas_rs_pool=auxbas_rs_pool,&
                                              para_env=para_env,&
                                              eps_mm_rspace=qmmm_env%eps_mm_rspace,&
                                              cube_info=ks_qmmm_env_loc%cube_info,&
                                              l_info=ks_qmmm_env_loc%l_info,&
                                              error=error)
       qs_env%ks_qmmm_env => ks_qmmm_env_loc       
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
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,error)
    Implicit None
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_electrostatic_on_grid', &
         routineP = moduleN//':'//routineN
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms

    CALL qmmm_electrostatic_on_grid_low(v_qmmm,v_qmmm%pw%cr3d,&
       mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,error)
  END SUBROUTINE qmmm_electrostatic_on_grid

!!****f* qmmm_electrostatic_methods/qmmm_electrostatic_on_grid_low [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_on_grid_low
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction on Grid Points: Low Level
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
  SUBROUTINE  qmmm_electrostatic_on_grid_low(v_qmmm,Varray,mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,error)
    Implicit None
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    REAL(dbl), DIMENSION(*) :: Varray
    INTEGER, INTENT(IN) :: num_mm_atoms
    ! Local Variables
    LOGICAL :: failure
    INTEGER :: handle
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER                   :: n1,n2,n3
    INTEGER :: i, j, k, Imm, IndMM,jj,my_k,my_j
    REAL(dbl) :: r, r0, qt, rNexp, rNexp1, r0Nexp, r0Nexp1
    REAL(dbl) :: dr1,dr2,dr3, xs1,xs2,xs3, rv1,rv2,rv3, rt1,rt2,rt3
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_electrostatic_on_grid', &
         routineP = moduleN//':'//routineN
    ! Statements
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)

    bo = v_qmmm%pw%pw_grid%bounds_local  
    n1 = v_qmmm%pw%pw_grid%npts (1)   
    n2 = v_qmmm%pw%pw_grid%npts (2)   
    n3 = v_qmmm%pw%pw_grid%npts (3)   
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
       qt=mm_charges(Imm)
       !$omp parallel do default (none) private( i, j, k, Imm, IndMM,jj,my_k,my_j,&
       !$omp                            r, r0, qt, rNexp, rNexp1, r0Nexp, r0Nexp1,&
       !$omp                   dr1,dr2,dr3, xs1,xs2,xs3, rv1,rv2,rv3, rt1,rt2,rt3)&
       !$omp                            shared(Varray,bo,n1,n2,n3)
       LoopOnGrid: Do k = bo(1,3), bo(2,3)
          my_k=(k+n3/2)
          xs3 = real(my_k,dbl)*dr3
          my_j=bo(1,2)+n2/2
          jj=my_k*n1*n2+my_j*n1
          xs2 = real(my_j,dbl)*dr2
          rv3 = rt3 - xs3
          Do j = bo(1,2), bo(2,2)
             xs1 = (bo(1,1)+n1/2)*dr1
             rv2 = rt2 - xs2
             Do i = bo(1,1), bo(2,1)
                jj=jj+1
                rv1 = rt1 - xs1
                r  = sqrt(rv1*rv1 + rv2*rv2 + rv3*rv3)
                rNexp = r**2
                rNexp = rNexp*rNexp
                rNexp1 = rNexp*r
                Varray(jj)=Varray(jj)+ &
                                        (-1.0_dbl) * & ! Sign of electrons
                                        qt* &
                                        (rNexp - r0Nexp)/(rNexp1 - r0Nexp1)
                xs1 = xs1 + dr1
             End Do
             xs2 = xs2 + dr2
          End Do
       End Do LoopOnGrid
    End Do MainLoopMM
    CALL timestop(0.0_dp,handle)

  End SUBROUTINE qmmm_electrostatic_on_grid_low

!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_with gaussian
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating the gaussian
!!     Electrostatic Potential
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
!!     06.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian(v_qmmm, mm_particles,&
       mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
       auxbas_rs_pool, cube_info, l_info, para_env, eps_mm_rspace, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    REAL(KIND=dbl), INTENT(IN) :: eps_mm_rspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), POINTER            :: cube_info
    TYPE(l_info_type)   , POINTER            :: l_info    
    ! Local Variables
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian', &
         routineP = moduleN//':'//routineN
    TYPE(realspace_grid_type), POINTER       :: tmp_grid
    ! Statements
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)
    !Statements
    NULLIFY(tmp_grid)
    CALL rs_pool_create_rs(auxbas_rs_pool,tmp_grid)
    CALL rs_grid_zero(tmp_grid)
    
    CALL qmmm_elec_with_gaussian_low( v_qmmm, tmp_grid, mm_particles,&
         mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms, cube_info,&
         l_info, para_env, eps_mm_rspace, error)

    CALL rs_pool_give_back_rs(auxbas_rs_pool,tmp_grid,error=error)

  END SUBROUTINE qmmm_elec_with_gaussian

!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian_low [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_with gaussian_low
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating the gaussian
!!     Electrostatic Potential - Low Level
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
!!     06.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_low( v_qmmm, tmp_grid, mm_particles,&
         mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
         cube_info, l_info, para_env, eps_mm_rspace, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(realspace_grid_type), POINTER           :: tmp_grid
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms
    TYPE(cube_info_type), POINTER         :: cube_info
    TYPE(l_info_type)   , POINTER         :: l_info
    REAL(KIND=dbl), INTENT(IN)               :: eps_mm_rspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_low', &
         routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    INTEGER, DIMENSION(:, :), POINTER        :: tasks
    LOGICAL :: failure
    INTEGER :: ncurr, stat, handle, ithread, j
    INTEGER :: dir, npme, Imm, IndMM
    INTEGER, DIMENSION(2) :: bo
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp) :: alpha
    !Statements
    NULLIFY(pab, tasks)
    CALL timeset(routineN,"I","",handle)
    ALLOCATE(pab(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ncurr = -1

    bo = get_limit ( num_mm_atoms, para_env%num_pe, para_env%mepos )
    npme = bo(2) - bo(1) + 1
    CALL reallocate ( tasks, 1, 2, 1, npme )
    ncurr = npme
    dir = tmp_grid%direction

!$OMP parallel do &
!$OMP default(none) &
!$OMP private(Imm,j,IndMM,ra) &
!$OMP shared(bo,tasks,dir,tmp_grid)

    DO Imm = bo(1), bo(2)
       j = Imm - bo(1) + 1
       tasks ( 1, j ) = Imm
       IF ( dir /= 0) THEN
          IndMM = mm_atom_index(Imm)
          ra(:) = mm_particles(IndMM)%r
          !tasks ( 2, j ) = FLOOR(ra(dir)/tmp_grid%dr(dir))
          !tasks ( 2, j ) = MODULO ( tasks ( 2, j ), tmp_grid%npts(dir) )
          !tasks ( 2, j ) = tasks ( 2, j ) + tmp_grid%lb(dir)
       END IF
    END DO
    
    CALL rs_get_my_tasks ( tmp_grid, tasks, npme )

    ithread=0
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,IndMM,Imm,ra,ithread) &
!$OMP shared(npme,tasks,mm_atom_index,mm_particles,pab,tmp_grid) &
!$OMP shared(cube_info,l_info,eps_mm_rspace)
    DO j = 1, npme
       Imm   = tasks(1,j)
       IndMM = mm_atom_index(Imm)
       ra(:) = mm_particles(IndMM)%r
       alpha     = 1.0_dbl / ( 2.0_dbl * mm_el_pot_radius(Imm)**2 )
       pab(1,1)  = - mm_charges(Imm) * SQRT(alpha/pi)**3
       WRITE(*,'(A,I5,A,I5,A,3F12.6,A,F12.6,A,F12.6)')&
                        ' Collocation Number :',j,&
                        ' MM Atom:',IndMM,&
                        ' x,y,z:',ra(1), ra(2), ra(3),&
                        ' alpha:',alpha,&
                        ' prefactor:',pab(1,1)
!$     ithread=omp_get_thread_num()
       CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dp,0,ra,&
            (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,-1.0_dp,pab,0,0,tmp_grid,&
            cube_info,l_info,eps_mm_rspace,ithread=ithread)
    END DO

    IF (ASSOCIATED(tasks)) THEN
       DEALLOCATE (tasks,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(pab)) THEN
       DEALLOCATE (pab, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL rs_pw_transfer(tmp_grid, v_qmmm%pw,"FORWARD")
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_elec_with_gaussian_low

End MODULE qmmm_electrostatic_methods
