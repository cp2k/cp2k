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
  USE global_types,                    ONLY: global_environment_type
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
                                             pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
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
  USE util,                            ONLY: get_limit
  USE qs_collocate_density,            ONLY: calculate_total_rho
  USE mm_collocate_potential,          ONLY: collocate_gf_rspace_NoPBC
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE pw_types,                        ONLY: pw_type,&
                                             pw_zero,&
                                             pw_prolongate_NoPBC
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_electrostatic_methods'

  PUBLIC :: qmmm_el_coupling
  PUBLIC :: qmmm_electrostatic_on_grid
  PUBLIC :: qmmm_elec_with_gaussian, qmmm_elec_with_gaussian_LR
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
  SUBROUTINE  qmmm_el_coupling(globenv,qs_env,qmmm_env,mm_particles,error)
    Implicit None
    ! Arguments
    TYPE(global_environment_type), &
      INTENT(in), TARGET,  OPTIONAL              :: globenv
    TYPE(qs_environment_type), INTENT(INOUT), &
         TARGET                                  :: qs_env   
    TYPE(qmmm_env_type), INTENT(in), &
         TARGET                                  :: qmmm_env  
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_el_coupling', &
         routineP = moduleN//':'//routineN    
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    LOGICAL :: failure
    INTEGER :: handle
    TYPE(pw_pool_p_type), POINTER,DIMENSION(:)     :: pw_pools
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=80)  :: FileName
    CHARACTER(LEN=80)  :: Title
    ! Statements
    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(ks_qmmm_env_loc, pw_pools, pw_env)
    CALL get_qs_env(qs_env=qs_env,&
                    pw_env=pw_env,&
                    para_env=para_env,&
                    error=error)

    CALL pw_env_get(pw_env=pw_env, pw_pools=pw_pools, error=error)

    ks_qmmm_env_loc => qs_env%ks_qmmm_env
    Write(6,'(T3,A)')"Information on the QM/MM Electrostatic Potential:"
    !
    ! Initializing vectors:
    !        Zeroing v_qmmm_rspace 
    CALL coeff_zero(coeff=ks_qmmm_env_loc%v_qmmm_rspace)
    !
    ! Electrostatic Interaction type...
    !
    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE("GRID")
       Write(6,'(T3,A)')"QM/MM Coupling computed on GRID."
       ! Compute the MM electrostatic potential on GRID
       CALL qmmm_electrostatic_on_grid(v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                       mm_particles=mm_particles,&
                                       mm_charges=qmmm_env%mm_atom_chrg,&
                                       mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                       mm_atom_index=qmmm_env%mm_atom_index,&
                                       num_mm_atoms=qmmm_env%num_mm_atoms,&
                                       error=error)
    CASE("GAUSS")
       Write(6,'(T3,A)')"QM/MM Coupling computed collocating the Gaussian Potential Functions."
       CALL    qmmm_elec_with_gaussian(v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                       mm_particles=mm_particles,&
                                       mm_charges=qmmm_env%mm_atom_chrg,&
                                       mm_el_pot_radius=qmmm_env%mm_el_pot_radius,&
                                       mm_atom_index=qmmm_env%mm_atom_index,&
                                       num_mm_atoms=qmmm_env%num_mm_atoms,&
                                       aug_pools=qmmm_env%aug_pools,&
                                       para_env=para_env,&
                                       eps_mm_rspace=qmmm_env%eps_mm_rspace,&
                                       cube_info=ks_qmmm_env_loc%cube_info,&
                                       pgfs=qmmm_env%pgfs,&
                                       pw_pools=pw_pools,&
                                       auxbas_grid=qmmm_env%gridlevel_info%auxbas_grid,&
                                       coarser_grid=qmmm_env%gridlevel_info%coarser_grid,&
                                       error=error)
    CASE("NONE")
       Write(6,'(T3,A)')"No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
    CASE DEFAULT
       Write(6,'(T3,A)')"Unknown Coupling..."
       CPPrecondition(.false.,cp_failure_level,routineP,error,failure)
    END SELECT

    qs_env%ks_qmmm_env => ks_qmmm_env_loc

    IF (globenv%print%v_qmmm_cube) THEN
       FileName = "MM_ELPOT_"//TRIM(qmmm_env%qmmm_coupl_type)//".cube"
       Title    = "QM/MM: MM ELECTROSTATIC POTENTIAL "
       IF (globenv%ionode) OPEN(UNIT=21,FILE=FileName)
       CALL rs_pw_to_cube(ks_qmmm_env_loc%v_qmmm_rspace%pw,21,globenv%ionode,Title)
       IF (globenv%ionode) CLOSE(21)
    END IF


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
    INTEGER :: handle

    CALL timeset(routineN,"I","",handle)
    CALL qmmm_electrostatic_on_grid_low(v_qmmm,v_qmmm%pw%cr3d,&
       mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,error)
    CALL timestop(0.0_dp,handle)
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
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_electrostatic_on_grid_low', &
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
       aug_pools, cube_info, para_env, eps_mm_rspace, pgfs, pw_pools, auxbas_grid,&
       coarser_grid, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(pw_pool_p_type), POINTER, DIMENSION(:) :: pw_pools, aug_pools
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, INTENT(IN) :: num_mm_atoms
    REAL(KIND=dbl), INTENT(IN) :: eps_mm_rspace
    INTEGER, INTENT(IN) :: auxbas_grid, coarser_grid
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                   :: pgfs
    ! Local Variables
    TYPE(pw_p_type), DIMENSION(:), pointer :: grids
    LOGICAL :: failure
    INTEGER :: ngrids,igrid, handle
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian', &
         routineP = moduleN//':'//routineN
    ! Statements
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(aug_pools),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pools),cp_failure_level,routineP,error,failure)
    !Statements
    CALL timeset(routineN,"I","",handle)
    ngrids=SIZE(pw_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw)
    END DO
 
    CALL qmmm_elec_with_gaussian_low( v_qmmm, grids, mm_particles,&
         mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms, cube_info,&
         para_env, eps_mm_rspace, pgfs, auxbas_grid, coarser_grid, pw_pools, aug_pools, &
         error)

    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)

    CALL timestop(0.0_dp,handle)

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
         cube_info, para_env, eps_mm_rspace, pgfs, auxbas_grid, coarser_grid,&
         pw_pools, aug_pools, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                :: error
    TYPE(pw_pool_p_type), POINTER, DIMENSION(:) :: pw_pools, aug_pools
    TYPE(coeff_type), INTENT(INOUT)              :: v_qmmm
    TYPE(pw_p_type), DIMENSION(:), POINTER       :: tmp_grid
    TYPE(particle_type), DIMENSION(:), POINTER   :: mm_particles
    REAL(dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    INTEGER, INTENT(IN) :: num_mm_atoms, auxbas_grid, coarser_grid
    TYPE(cube_info_type), DIMENSION(:), POINTER         :: cube_info
    REAL(KIND=dbl), INTENT(IN)            :: eps_mm_rspace
    TYPE(cp_para_env_type), POINTER       :: para_env
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                   :: pgfs
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_low', &
         routineP = moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat, handle, ithread, ilevel, igrid
    INTEGER :: Imm, IndMM, LIndMM, IRadTyp, Icg, IGauss
    INTEGER, DIMENSION(3) :: ndim, lb, ub
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp) :: alpha, height, ConsFac, W
    TYPE(qmmm_gaussian_type), POINTER :: pgf
    REAL(KIND=dbl), DIMENSION(:), POINTER    :: xdat, ydat, zdat
    !Statements
    NULLIFY(pgf, xdat, ydat, zdat)
    CALL timeset(routineN,"I","",handle)

    ithread            = 0
    Icg                = 0
    ndim               = 0
    DO igrid = auxbas_grid, coarser_grid
       ndim = MAX(tmp_grid(igrid)%pw%pw_grid%npts,ndim)
    END DO
    ALLOCATE (xdat(ndim(1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(ndim(2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(ndim(3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Collocating the Constant Term (Just for Energy Consistency)       
    ConsFac=0._dbl
    DO IRadTyp = 1, SIZE(pgfs) 
       pgf => pgfs(IRadTyp)%pgf
       DO Imm  = 1, SIZE(pgf%mm_atom_index)
          LIndMM    =   pgf%mm_atom_index(Imm)
          ConsFac   = - pgf%A0 * mm_charges(LIndMM) + ConsFac
       END DO
    END DO
    tmp_grid(coarser_grid)%pw%cr3d = ConsFac

!
    Radius: DO IRadTyp = 1, SIZE(pgfs) 
       pgf => pgfs(IRadTyp)%pgf
       !WRITE(*,'(A,F12.6)')' Collocating MM Potential with Radius: ', pgf%Elp_Radius
       Gaussian: DO IGauss = 1, pgf%Number_of_Gaussians
          alpha     = 1.0_dbl / pgf%Gk(IGauss)
          alpha     = alpha * alpha
          height    = pgf%Ak(IGauss)
          ilevel    = pgf%grid_level(IGauss)
          IF (ilevel.EQ.coarser_grid) CYCLE Gaussian
          !WRITE(*,'(5X,A,I5)')' Collocating Gaussian Number:',IGauss
          Atoms: DO Imm = 1, SIZE(pgf%mm_atom_index)
             Icg       =   Icg + 1
             LIndMM    =   pgf%mm_atom_index(Imm)
             IndMM     =   mm_atom_index(LIndMM)
             ra(:)     =   mm_particles(IndMM)%r
             W         =   mm_charges(LIndMM) * height
             !WRITE(*,'(10X,I5,A,I5,A,3F12.6,A,F12.6,A,F12.6)')&
             !     Icg,' MM Atom:',IndMM,&
             !     ' x,y,z:',ra(1), ra(2), ra(3),&
             !     ' alpha:',alpha,&
             !     ' prefactor:',W
             
             CALL collocate_gf_rspace_NoPBC(zetp=alpha,&
                                            rp=ra,&
                                            scale=-1.0_dp,&
                                            W=W,&
                                            pwgrid=tmp_grid(ilevel)%pw,&
                                            cube_info=cube_info(ilevel),&
                                            eps_mm_rspace=eps_mm_rspace,&
                                            ithread=ithread, &
                                            xdat=xdat,&
                                            ydat=ydat,&
                                            zdat=zdat)                          
          END DO Atoms
       END DO Gaussian
    END DO Radius

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

    ! Long Range Part of the QM/MM Potential with Gaussians
    CALL qmmm_elec_with_gaussian_LR    (pgfs=pgfs,&
                                        grid=tmp_grid(coarser_grid)%pw,&
                                        mm_charges=mm_charges,&
                                        mm_atom_index=mm_atom_index,&
                                        mm_particles=mm_particles,&
                                        mm_el_pot_radius=mm_el_pot_radius,&
                                        para_env=para_env,&
                                        coarser_grid_level=coarser_grid,&
                                        error=error)

    ! RealSpace Interpolation

    IF (.true.) THEN
       ! Linear Interpolator
       DO Ilevel = coarser_grid, auxbas_grid+1, -1
          CALL pw_prolongate_NoPBC(tmp_grid(Ilevel  )%pw,&
                                   tmp_grid(Ilevel-1)%pw,&
                                   error)
       END DO
    ELSE
       ! Spline Interpolator - Still to be completed... 
       ! Fawzi?!?! Yuuuuuhuuuuuu?!?! 
       !
       !DO Ilevel = coarser_grid, auxbas_grid+1, -1
       !   CALL pw_inject_s3(tmp_grid(Ilevel  )%pw,&
       !                     tmp_grid(Ilevel-1)%pw,&
       !                     pw_pools(Ilevel-1)%pool,&
       !                     aug_pools(Ilevel-1)%pool,&
       !                     error)
       !END DO       
    ENDIF
    lb = v_qmmm%pw%pw_grid%bounds(1,:)
    ub = v_qmmm%pw%pw_grid%bounds(2,:)

    v_qmmm%pw%cr3d = tmp_grid(auxbas_grid)%pw%cr3d (lb(1):ub(1),&
                                                    lb(2):ub(2),&
                                                    lb(3):ub(3) )

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_elec_with_gaussian_low

!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian_LR [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_with gaussian_LR
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating 
!!     (1/R - Sum_NG Gaussians) on the coarser grid level.
!!     Long Range QM/MM Electrostatic Potential with Gaussian - Low Level
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
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_LR(pgfs, grid, mm_charges,&
       mm_atom_index, mm_particles, mm_el_pot_radius, para_env,&
       coarser_grid_level, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                  :: error
    TYPE(pw_type), POINTER                         :: grid
    TYPE(particle_type), DIMENSION(:), POINTER     :: mm_particles
    REAL(KIND=dbl), DIMENSION(:), POINTER :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER   :: mm_atom_index
    TYPE(cp_para_env_type), POINTER       :: para_env
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                   :: pgfs
    INTEGER, INTENT(IN) :: coarser_grid_level
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_LR', &
         routineP = moduleN//':'//routineN    
    TYPE(qmmm_gaussian_type), POINTER :: pgf
    INTEGER :: handle, IRadTyp, IGauss, Imm
    INTEGER :: IndMM, LIndMM, Ilevel
    INTEGER :: i, j, k, my_k, my_j, n1, n2, n3
    INTEGER, DIMENSION(2,3)        :: bo
    REAL(KIND=dbl), DIMENSION(3)   :: ra
    REAL(KIND=dbl)               :: W, height, alpha, Gterm
    REAL(KIND=dbl)               :: qt, rt1, rt2, rt3, r, r0, r2, r0Nexp, r0Nexp1
    REAL(KIND=dbl)               :: xs1, rv1, xs2, rv2, xs3, rv3, rNexp, rNexp1, dr1, dr2, dr3
    ! Statements
    CALL timeset(routineN,"I","",handle)

    n1  = grid%pw_grid%npts(1)
    n2  = grid%pw_grid%npts(2)
    n3  = grid%pw_grid%npts(3)
    dr1 = grid%pw_grid%dr(1) 
    dr2 = grid%pw_grid%dr(2) 
    dr3 = grid%pw_grid%dr(3) 
    bo(1,:) = grid%pw_grid%bounds(1,:)
    bo(2,:) = grid%pw_grid%bounds(2,:)

    Radius: DO IRadTyp = 1, SIZE(pgfs) 
       pgf => pgfs(IRadTyp)%pgf
       r0      = pgf%Elp_Radius
       r0Nexp  = r0**2
       r0Nexp  = r0Nexp*r0Nexp
       r0Nexp1 = r0*r0Nexp
       Atoms: DO Imm = 1, SIZE(pgf%mm_atom_index)
          LIndMM    =   pgf%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   mm_particles(IndMM)%r
          qt        =   mm_charges(LIndMM)
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k = (k+n3/2)
             xs3  = REAL(my_k,dbl)*dr3
             my_j = bo(1,2)+n2/2
             xs2 = REAL(my_j,dbl)*dr2
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                xs1 = (bo(1,1)+n1/2)*dr1
                rv2 = rt2 - xs2
                DO i =  bo(1,1), bo(2,1)
                   rv1 = rt1 - xs1
                   r2  = rv1*rv1 + rv2*rv2 + rv3*rv3
                   r   = SQRT(r2)
                   rNexp  = r**2
                   rNexp  = rNexp*rNexp
                   rNexp1 = rNexp*r
                   Gterm  = 0.0_dbl
                   Gaussian: DO IGauss = 1, pgf%Number_of_Gaussians
                      ilevel    = pgf%grid_level(IGauss)
                      IF (ilevel.EQ.coarser_grid_level) CYCLE Gaussian
                      alpha  = 1.0_dbl / pgf%Gk(IGauss)
                      alpha  = alpha * alpha
                      height = pgf%Ak(IGauss)
                      W      = qt * height
                      Gterm  = Gterm + W * EXP(-alpha*r2)
                   END DO Gaussian
                   Gterm = Gterm + pgf%A0 * qt
                   grid%cr3d(i,j,k) = grid%cr3d(i,j,k)    + &
                                  (-1.0_dbl) * qt * &
                                  (rNexp - r0Nexp)/(rNexp1 - r0Nexp1) + Gterm

                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
       END DO Atoms
    END DO Radius

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_elec_with_gaussian_LR


End MODULE qmmm_electrostatic_methods
