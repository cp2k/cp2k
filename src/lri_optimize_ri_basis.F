!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Optimizes exponents and contraction coefficients of the lri auxiliary
!>        basis sets using the UOBYQA minimizer 
!>        lri : local resolution of the identity
!> \par History
!>      created Dorothea Golze [05.2014]
!> \authors Dorothea Golze  
! *****************************************************************************
MODULE lri_optimize_ri_basis

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type,&
                                             init_orb_basis_set
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_lri_opt_all,&
                                             do_lri_opt_coeff,&
                                             do_lri_opt_exps
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE lri_environment_init,            ONLY: lri_basis_init
  USE lri_environment_methods,         ONLY: calculate_avec,&
                                             calculate_lri_integrals,&
                                             calculate_lri_overlap_aabb
  USE lri_environment_types,           ONLY: lri_density_type,&
                                             lri_environment_type,&
                                             lri_int_rho_type,&
                                             lri_int_type,&
                                             lri_list_type,&
                                             lri_rhoab_type
  USE lri_optimize_ri_basis_types,     ONLY: create_lri_opt,&
                                             deallocate_lri_opt,&
                                             get_original_gcc,&
                                             lri_opt_type,&
                                             orthonormalize_gcc
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE powell,                          ONLY: opt_state_type,&
                                             powell_optimize
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_optimize_ri_basis'

  PUBLIC :: optimize_lri_basis

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief optimizes the lri basis set
!> \param qs_env qs environment
! *****************************************************************************
  SUBROUTINE optimize_lri_basis(qs_env)

    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_lri_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iunit, nkind
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_opt_type), POINTER              :: lri_opt
    TYPE(opt_state_type)                     :: opt_state
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                lri_optbas_section

    NULLIFY(atomic_kind_set, dft_section, lri_density, lri_env, &
           lri_opt,lri_optbas_section,  rho_struct)
    NULLIFY(input, logger, para_env)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,input=input,&
                    lri_env=lri_env,lri_density=lri_density,nkind=nkind,&
                    para_env=para_env,rho=rho_struct)

    ! density matrix
    CALL qs_rho_get(rho_struct, rho_ao=pmatrix)

    logger => cp_get_default_logger()
    dft_section => section_vals_get_subs_vals(input,"DFT")
    lri_optbas_section => section_vals_get_subs_vals(input,&
                          "DFT%QS%OPTIMIZE_LRI_BASIS")
    iunit=cp_print_key_unit_nr(logger,input,"PRINT%PROGRAM_RUN_INFO",&
                               extension=".opt")

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")')
    ENDIF

    ! *** initialization
    CALL create_lri_opt(lri_opt)
    CALL init_optimization(lri_env,lri_opt,lri_optbas_section,&
                           opt_state,lri_opt%x,lri_opt%zet_init,nkind,iunit)

    CALL calculate_lri_overlap_aabb(lri_env,qs_env)

    ! *** ======================= START optimization ===================== 
    opt_state%state = 0
    DO
      IF ( opt_state%state == 2 ) THEN
        CALL calc_lri_integrals_get_objective(lri_env,lri_density,qs_env,&
                                              lri_opt,opt_state,pmatrix,para_env,&
                                              nkind)
      ENDIF

      IF ( opt_state%state == -1 ) EXIT

      CALL powell_optimize (opt_state%nvar, lri_opt%x, opt_state)
      CALL update_exponents(lri_env,lri_opt,lri_opt%x,lri_opt%zet_init,nkind)
      CALL print_optimization_update(opt_state,lri_opt,iunit)
    ENDDO
    ! *** ======================= END optimization =======================

    ! *** get final optimized parameters
    opt_state%state = 8
    CALL powell_optimize (opt_state%nvar, lri_opt%x, opt_state)
    CALL update_exponents(lri_env,lri_opt,lri_opt%x,lri_opt%zet_init,nkind)

    CALL write_optimized_lri_basis(lri_env,dft_section,nkind,lri_opt,&
         atomic_kind_set)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') opt_state%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,F20.10)') opt_state%fopt
      WRITE(iunit,'(/," Printed optimized lri basis set to file")')
    ENDIF

    CALL cp_print_key_finished_output(iunit,logger,input,&
                    "PRINT%PROGRAM_RUN_INFO")

    CALL deallocate_lri_opt(lri_opt)

  END SUBROUTINE optimize_lri_basis

! *****************************************************************************
!> \brief initialize optimization parameter
!> \param lri_env lri environment
!> \param lri_opt optimization environment
!> \param lri_optbas_section ...
!> \param opt_state state of the optimizer
!> \param x parameters to be optimized, i.e. exponents and contraction coeffs
!>        of the lri basis set 
!> \param zet_init initial values of the exponents
!> \param nkind number of atom kinds
!> \param iunit output unit
! *****************************************************************************
  SUBROUTINE init_optimization(lri_env,lri_opt,lri_optbas_section,opt_state,&
                               x,zet_init,nkind,iunit)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_opt_type), POINTER              :: lri_opt
    TYPE(section_vals_type), POINTER         :: lri_optbas_section
    TYPE(opt_state_type)                     :: opt_state
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x, zet_init
    INTEGER, INTENT(IN)                      :: nkind, iunit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_optimization', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, iset, ishell, n, nset, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: npgf, nshell
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc_orig
    TYPE(gto_basis_set_type), POINTER        :: fbas

    NULLIFY(fbas, gcc_orig, npgf, nshell, zet)
    failure = .FALSE.

    ALLOCATE(lri_opt%ri_gcc_orig(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    ! *** get parameters 
    CALL get_optimization_parameter(lri_opt,lri_optbas_section,&
                                    opt_state)

    opt_state%nvar =0
    opt_state%nf = 0
    opt_state%iprint = 1
    opt_state%unit  = iunit

    ! *** init exponents
    IF(lri_opt%opt_exps) THEN
      n = 0
      DO ikind=1,nkind
        fbas  => lri_env%ri_basis(ikind)%gto_basis_set
        CALL get_gto_basis_set(gto_basis_set=fbas,&
             npgf=npgf,nset=nset,zet=zet)
        DO iset =1,nset
          IF(lri_opt%use_geometric_seq.AND.npgf(iset)>2) THEN
            opt_state%nvar = opt_state%nvar + 2
            CALL reallocate(x,1,opt_state%nvar)
            x(n+1) = MAXVAL(zet(1:npgf(iset),iset))
            x(n+2) = MINVAL(zet(1:npgf(iset),iset))
            n = n + 2
          ELSE
            opt_state%nvar = opt_state%nvar + npgf(iset)
            CALL reallocate(x,1,opt_state%nvar)
            x(n+1:n+npgf(iset)) = zet(1:npgf(iset),iset)
            n = n + npgf(iset)
          ENDIF
          lri_opt%nexp=lri_opt%nexp + npgf(iset)
        ENDDO
      ENDDO
     
      ! *** constraints on exponents
      IF(lri_opt%use_constraints) THEN
        ALLOCATE(zet_init(SIZE(x)),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
        zet_init(:)=x 
      ELSE
        x(:)=SQRT(x)
      ENDIF
    ENDIF

    ! *** get the original gcc without normalization factor
    DO ikind=1,nkind
      fbas => lri_env%ri_basis(ikind)%gto_basis_set
      CALL get_original_gcc(lri_opt%ri_gcc_orig(ikind)%gcc_orig,fbas,&
                            lri_opt)
    ENDDO

    ! *** init coefficients
    IF(lri_opt%opt_coeffs) THEN
      DO ikind=1,nkind
        fbas => lri_env%ri_basis(ikind)%gto_basis_set
        gcc_orig => lri_opt%ri_gcc_orig(ikind)%gcc_orig
        CALL get_gto_basis_set(gto_basis_set=fbas,&
             npgf=npgf,nset=nset,nshell=nshell,zet=zet)
        ! *** Gram Schmidt orthonormalization
        CALL orthonormalize_gcc(gcc_orig,fbas,lri_opt)
        n=opt_state%nvar
        DO iset=1,nset
          DO ishell=1,nshell(iset)
            opt_state%nvar = opt_state%nvar + npgf(iset)
            CALL reallocate(x,1,opt_state%nvar)
            x(n+1:n+npgf(iset))=gcc_orig(1:npgf(iset),ishell,iset)
            lri_opt%ncoeff = lri_opt%ncoeff + npgf(iset)
            n = n + npgf(iset) 
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF(iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Accuracy",T69,ES12.5)') opt_state%rhoend
      WRITE(iunit,'(" POWELL| Initial step size",T69,ES12.5)') opt_state%rhobeg
      WRITE(iunit,'(" POWELL| Maximum number of evaluations",T71,I10)')&
                                                               opt_state%maxfun
      WRITE(iunit,'(" POWELL| Total number of parameters",T71,I10)')&
                                                               opt_state%nvar
    ENDIF

  END SUBROUTINE init_optimization

! *****************************************************************************
!> \brief read input for optimization
!> \param lri_opt optimization environmnet
!> \param lri_optbas_section ...
!> \param opt_state state of the optimizer
! *****************************************************************************
  SUBROUTINE get_optimization_parameter(lri_opt,lri_optbas_section,&
                                       opt_state)

    TYPE(lri_opt_type), POINTER              :: lri_opt
    TYPE(section_vals_type), POINTER         :: lri_optbas_section
    TYPE(opt_state_type)                     :: opt_state

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_optimization_parameter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: degree_freedom
    TYPE(section_vals_type), POINTER         :: constrain_exp_section

    NULLIFY(constrain_exp_section)

    ! *** parameter for POWELL optimizer
    CALL section_vals_val_get(lri_optbas_section,"ACCURACY",&
                              r_val=opt_state%rhoend)
    CALL section_vals_val_get(lri_optbas_section,"STEP_SIZE",&
                              r_val=opt_state%rhobeg)
    CALL section_vals_val_get(lri_optbas_section,"MAX_FUN",&
                              i_val=opt_state%maxfun)

    ! *** parameters which are optimized, i.e. exps or coeff or both
    CALL section_vals_val_get(lri_optbas_section,"DEGREES_OF_FREEDOM",&
                              i_val=degree_freedom)

    SELECT CASE(degree_freedom)
    CASE(do_lri_opt_all)
       lri_opt%opt_coeffs=.TRUE.
       lri_opt%opt_exps=.TRUE.
    CASE(do_lri_opt_coeff)
       lri_opt%opt_coeffs=.TRUE.
    CASE(do_lri_opt_exps)
       lri_opt%opt_exps=.TRUE.
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,&
                 routineP,"No initialization available?????",&
                 only_ionode=.TRUE.)
    END SELECT

    ! *** restraint
    CALL section_vals_val_get(lri_optbas_section,"USE_CONDITION_NUMBER",&
                              l_val=lri_opt%use_condition_number)
    CALL section_vals_val_get(lri_optbas_section,"CONDITION_WEIGHT",&
                              r_val=lri_opt%cond_weight)
    CALL section_vals_val_get(lri_optbas_section,"GEOMETRIC_SEQUENCE",&
                              l_val=lri_opt%use_geometric_seq)

    ! *** get constraint info
    constrain_exp_section  => section_vals_get_subs_vals(lri_optbas_section,&
                              "CONSTRAIN_EXPONENTS")
    CALL section_vals_get(constrain_exp_section,explicit=lri_opt%use_constraints)

    IF(lri_opt%use_constraints) THEN
      CALL section_vals_val_get(constrain_exp_section,"SCALE",&
                                r_val=lri_opt%scale_exp)
      CALL section_vals_val_get(constrain_exp_section,"FERMI_EXP",&
                                r_val=lri_opt%fermi_exp)
    ENDIF

  END SUBROUTINE get_optimization_parameter

! *****************************************************************************
!> \brief update exponents after optimization step
!> \param lri_env lri environment
!> \param lri_opt optimization environment
!> \param x optimization parameters
!> \param zet_init initial values of the exponents
!> \param nkind number of atomic kinds
! *****************************************************************************
  SUBROUTINE update_exponents(lri_env,lri_opt,x,zet_init,nkind)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_opt_type), POINTER              :: lri_opt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x, zet_init
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_exponents', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, iset, ishell, n, nset, &
                                                nvar_exp, stat
    INTEGER, DIMENSION(:), POINTER           :: npgf, nshell
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zet_max, zet_min
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet, zet_trans
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc_orig
    TYPE(gto_basis_set_type), POINTER        :: fbas

    NULLIFY(fbas,gcc_orig,npgf,nshell,zet_trans, zet)
    failure = .FALSE.

    ! nvar_exp: number of exponents that are variables
    nvar_exp= SIZE(x) - lri_opt%ncoeff
    ALLOCATE(zet_trans(nvar_exp),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
   
    ! *** update exponents 
    IF(lri_opt%opt_exps) THEN
      IF (lri_opt%use_constraints) THEN
        zet => x(1:nvar_exp)
        CALL transfer_exp(lri_opt,zet,zet_init,zet_trans,nvar_exp)
      ELSE
        zet_trans(:) = x(1:nvar_exp)**2.0_dp
      ENDIF
      n=0
      DO ikind =1,nkind
        fbas  => lri_env%ri_basis(ikind)%gto_basis_set
        CALL get_gto_basis_set(gto_basis_set=fbas,npgf=npgf,nset=nset)
        DO iset=1,nset
          IF(lri_opt%use_geometric_seq.AND.npgf(iset)>2) THEN
            zet_max = MAXVAL(zet_trans(n+1:n+2))
            zet_min = MINVAL(zet_trans(n+1:n+2))
            zet => fbas%zet(1:npgf(iset),iset)
            CALL geometric_progression(zet,zet_max,zet_min,npgf(iset))
            n = n + 2
          ELSE
            fbas%zet(1:npgf(iset),iset) = zet_trans(n+1:n+npgf(iset))
            n = n + npgf(iset)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
   
    ! *** update coefficients
    IF(lri_opt%opt_coeffs) THEN
      n=nvar_exp
      DO ikind=1,nkind
        fbas  => lri_env%ri_basis(ikind)%gto_basis_set
        gcc_orig => lri_opt%ri_gcc_orig(ikind)%gcc_orig
        CALL get_gto_basis_set(gto_basis_set=fbas,&
             nshell=nshell,npgf=npgf,nset=nset)
        DO iset=1,nset
          DO ishell=1,nshell(iset)
            gcc_orig(1:npgf(iset),ishell,iset) = x(n+1:n+npgf(iset))
            n = n + npgf(iset)
          ENDDO
        ENDDO
        ! *** Gram Schmidt orthonormalization
        CALL orthonormalize_gcc(gcc_orig,fbas,lri_opt)
      ENDDO
    ENDIF

    DEALLOCATE(zet_trans,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
  END SUBROUTINE update_exponents

! *****************************************************************************
!> \brief employ Fermi constraint, transfer exponents
!> \param lri_opt optimization environment
!> \param zet untransferred exponents
!> \param zet_init intial value of the eponents
!> \param zet_trans transferred exponents
!> \param nvar number of optimized exponents 
! *****************************************************************************
  SUBROUTINE transfer_exp(lri_opt,zet,zet_init,zet_trans,nvar)

    TYPE(lri_opt_type), POINTER              :: lri_opt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet, zet_init, zet_trans
    INTEGER, INTENT(IN)                      :: nvar

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transfer_exp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet_max, zet_min

    failure = .FALSE.
    ALLOCATE(zet_max(nvar),zet_min(nvar),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    
    zet_min(:) = zet_init(:)*(1.0_dp-lri_opt%scale_exp)
    zet_max(:) = zet_init(:)*(1.0_dp+lri_opt%scale_exp)
   
    a=lri_opt%fermi_exp
   
    zet_trans= zet_min + (zet_max-zet_min)/(1+EXP(-a*(zet-zet_init)))

    DEALLOCATE(zet_max,zet_min,STAT=stat) 
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)


  END SUBROUTINE transfer_exp

! *****************************************************************************
!> \brief complete geometric sequence
!> \param zet all exponents of the set
!> \param zet_max maximal exponent of the set
!> \param zet_min minimal exponent of the set
!> \param nexp number of exponents of the set
! *****************************************************************************
  SUBROUTINE geometric_progression(zet,zet_max,zet_min,nexp)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet
    REAL(KIND=dp), INTENT(IN)                :: zet_max, zet_min
    INTEGER, INTENT(IN)                      :: nexp

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geometric_progression', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    REAL(KIND=dp)                            :: q

    n=nexp-1

    q=(zet_min/zet_max)**(1._dp/REAL(n,dp))

    DO i=1,nexp
     zet(i)=zet_max*q**(i-1)
    ENDDO

  END SUBROUTINE geometric_progression

! *****************************************************************************
!> \brief calculates the lri integrals and coefficients with the new exponents
!>        of the lri basis sets and calculates the objective function
!> \param lri_env lri environment
!> \param lri_density ...
!> \param qs_env ...
!> \param lri_opt optimization enviroment
!> \param opt_state state of the optimizer
!> \param pmatrix density matrix
!> \param para_env ...
!> \param nkind number of atomic kinds
! *****************************************************************************
  SUBROUTINE calc_lri_integrals_get_objective(lri_env,lri_density,qs_env,&
                                              lri_opt,opt_state,pmatrix,para_env,&
                                              nkind)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lri_opt_type), POINTER              :: lri_opt
    TYPE(opt_state_type)                     :: opt_state
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calc_lri_integrals_get_objective', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nset
    INTEGER, DIMENSION(:), POINTER           :: npgf
    TYPE(gto_basis_set_type), POINTER        :: fbas

    NULLIFY(fbas,npgf)

    !*** build new transformation matrices sphi with new exponents
    DO ikind =1,nkind
      fbas  => lri_env%ri_basis(ikind)%gto_basis_set
      CALL get_gto_basis_set(gto_basis_set=fbas,npgf=npgf,nset=nset)
      !build new sphi
      fbas%gcc = lri_opt%ri_gcc_orig(ikind)%gcc_orig
      CALL init_orb_basis_set(fbas)
    ENDDO
    CALL lri_basis_init(lri_env)
    CALL calculate_lri_integrals(lri_env,qs_env,calculate_forces=.FALSE.)
    CALL calculate_avec(lri_env,lri_density,qs_env,pmatrix)
    IF(lri_opt%use_condition_number) THEN 
      CALL get_condition_number_of_overlap(lri_env)
    ENDIF
    CALL calculate_objective(lri_env,lri_density,lri_opt,pmatrix,para_env,&
                             opt_state%f)


  END SUBROUTINE calc_lri_integrals_get_objective

! *****************************************************************************
!> \brief calculates the objective function defined as integral of the square
!>        of rhoexact - rhofit, i.e. integral[(rhoexact-rhofit)**2]
!>        rhoexact is the exact pair density and rhofit the lri pair density
!> \param lri_env lri environment
!> \param lri_density ...
!> \param lri_opt optimization environment
!> \param pmatrix density matrix
!> \param para_env ...
!> \param fobj objective function
! *****************************************************************************
  SUBROUTINE calculate_objective(lri_env,lri_density,lri_opt,pmatrix,para_env,&
                                 fobj)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_opt_type), POINTER              :: lri_opt
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(OUT)               :: fobj

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_objective', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iac, iatom, ikind, ilist, isgfa, ispin, jatom, jkind, &
      jneighbor, jsgfa, ksgfb, lsgfb, nba, nbb, nfa, nfb, nkind, nlist, nn, &
      nneighbor, nspin
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp)                            :: obj_ab, rhoexact_sq, &
                                                rhofit_sq, rhomix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pbij
    TYPE(cp_dbcsr_type), POINTER             :: pmat
    TYPE(lri_int_rho_type), POINTER          :: lriir
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_list_type), POINTER             :: lri_rho
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(lrii, lriir, lri_rho, lrho, nl_iterator, pmat, soo_list)

    IF (ASSOCIATED(lri_env%soo_list)) THEN
      soo_list => lri_env%soo_list

      nkind = lri_env%lri_ints%nkind
      nspin = SIZE(pmatrix)
      fobj = 0._dp
      lri_opt%rho_diff = 0._dp
       

      DO ispin = 1, nspin

        pmat => pmatrix(ispin)%matrix
        lri_rho => lri_density%lri_rhos(ispin)%lri_list

        CALL neighbor_list_iterator_create(nl_iterator,soo_list)
        DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,&
               jatom=jatom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor)

          iac = ikind + nkind*(jkind - 1)

          IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

          lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
          lriir => lri_env%lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)
          lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
          nfa = lrii%nfa
          nfb = lrii%nfb
          nba = lrii%nba
          nbb = lrii%nbb
          nn = nfa +  nfb
 
          rhoexact_sq = 0._dp 
          rhomix      = 0._dp 
          rhofit_sq   = 0._dp
          obj_ab      = 0._dp

          NULLIFY (pbij)
          IF (iatom <= jatom) THEN
            CALL cp_dbcsr_get_block_p(matrix=pmat,row=iatom,col=jatom,block=pbij,found=found)
            trans = .FALSE.
          ELSE
            CALL cp_dbcsr_get_block_p(matrix=pmat,row=jatom,col=iatom,block=pbij,found=found)
            trans = .TRUE.
          END IF
          CPPostcondition(found,cp_failure_level,routineP,failure)
 
          ! *** calculate integral of the square of exact density rhoexact_sq
          IF (trans) THEN  
            DO isgfa = 1, nba
              DO jsgfa = 1, nba
                DO ksgfb = 1, nbb
                  DO lsgfb = 1, nbb
                    rhoexact_sq = rhoexact_sq + pbij(ksgfb,isgfa)*pbij(lsgfb,jsgfa)&
                                               *lriir%soaabb(isgfa,jsgfa,ksgfb,lsgfb) 
                  END DO
                END DO
              ENDDO
            END DO  
          ELSE
            DO isgfa = 1, nba
              DO jsgfa = 1, nba
                DO ksgfb = 1, nbb
                  DO lsgfb = 1, nbb
                    rhoexact_sq = rhoexact_sq + pbij(isgfa,ksgfb)*pbij(jsgfa,lsgfb)&
                                               *lriir%soaabb(isgfa,jsgfa,ksgfb,lsgfb) 
                  END DO
                END DO
              ENDDO
            END DO  
          ENDIF

        
          ! *** calculate integral of the square of the fitted density rhofit_sq 
          DO isgfa=1,nfa
            DO jsgfa=1,nfa
              rhofit_sq = rhofit_sq + lrho%avec(isgfa)*lrho%avec(jsgfa)&
                                      *lri_env%bas_prop(ikind)%ri_ovlp(isgfa,jsgfa)
            ENDDO
          ENDDO
          IF(iatom /= jatom) THEN
            DO ksgfb=1,nfb
              DO lsgfb=1,nfb
                rhofit_sq = rhofit_sq + lrho%avec(nfa+ksgfb)*lrho%avec(nfa+lsgfb)&
                                       *lri_env%bas_prop(jkind)%ri_ovlp(ksgfb,lsgfb)
              ENDDO
            ENDDO
            DO isgfa=1,nfa
              DO ksgfb=1,nfb
                rhofit_sq= rhofit_sq + 2._dp * lrho%avec(isgfa)*lrho%avec(nfa+ksgfb)&
                                             * lrii%sab(isgfa,ksgfb) 
              ENDDO
            ENDDO
          ENDIF

          ! *** and integral of the product of exact and fitted density rhomix
          IF(iatom == jatom) THEN
            rhomix = SUM(lrho%avec(1:nfa) * lrho%tvec(1:nfa))
          ELSE
            rhomix = SUM(lrho%avec(1:nn) * lrho%tvec(1:nn))
          ENDIF

          ! *** calculate contribution to the objective function for pair ab
          ! *** taking density matrix symmetry in account, double-count for off-diagonal blocks
          IF(iatom == jatom) THEN
            obj_ab = rhoexact_sq - 2._dp*rhomix + rhofit_sq
          ELSE
            obj_ab = 2.0_dp*(rhoexact_sq - 2._dp*rhomix + rhofit_sq)
          ENDIF

          IF(lri_opt%use_condition_number) THEN
            fobj = fobj + obj_ab + lri_opt%cond_weight*LOG(lrii%cond_num)
            lri_opt%rho_diff = lri_opt%rho_diff + obj_ab
          ELSE
            fobj = fobj + obj_ab
          ENDIF

        ENDDO

        CALL neighbor_list_iterator_release(nl_iterator)

      ENDDO
      CALL mp_sum(fobj,para_env%group) 

    ENDIF 

    CALL timestop(handle)

  END SUBROUTINE calculate_objective

! *****************************************************************************
!> \brief get condition number of overlap matrix
!> \param lri_env lri environment
! *****************************************************************************
  SUBROUTINE get_condition_number_of_overlap(lri_env)

    TYPE(lri_environment_type), POINTER      :: lri_env

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'get_condition_number_of_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iac, iatom, ikind, ilist, info, jatom, jkind, &
      jneighbor, lwork, nfa, nfb, nkind, nlist, nn, nneighbor, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: diag, off_diag, tau
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: smat
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(lrii, nl_iterator, smat, soo_list)

    soo_list => lri_env%soo_list

    nkind = lri_env%lri_ints%nkind
    CALL neighbor_list_iterator_create(nl_iterator,soo_list)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
      CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,&
           jatom=jatom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor)

       iac = ikind + nkind*(jkind - 1)
       IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE
       lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

       nfa=lrii%nfa
       nfb=lrii%nfb
       nn = nfa + nfb

       ! build the overlap matrix
       IF(iatom == jatom) THEN
        ALLOCATE(smat(nfa,nfa),STAT=stat) 
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
       ELSE
        ALLOCATE(smat(nn,nn),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,failure)
       ENDIF
       smat(1:nfa,1:nfa) = lri_env%bas_prop(ikind)%ri_ovlp(1:nfa,1:nfa)
       IF(iatom /= jatom) THEN
         nn = nfa+nfb
         smat(1:nfa,nfa+1:nn) = lrii%sab(1:nfa,1:nfb)
         smat(nfa+1:nn,1:nfa) = TRANSPOSE(lrii%sab(1:nfa,1:nfb))
         smat(nfa+1:nn,nfa+1:nn) = lri_env%bas_prop(jkind)%ri_ovlp(1:nfb,1:nfb)
       ENDIF

       IF(iatom==jatom) nn=nfa
       ALLOCATE(diag(nn),off_diag(nn-1),tau(nn-1),work(1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,failure)
       diag=0.0_dp
       off_diag=0.0_dp
       tau=0.0_dp
       work=0.0_dp
       lwork = -1
       ! get lwork
       CALL DSYTRD('U',nn,smat,nn,diag,off_diag,tau,work,lwork,info)
       lwork=INT(work(1))
       CALL reallocate(work,1,lwork)
       ! get the eigenvalues
       CALL DSYTRD('U',nn,smat,nn,diag,off_diag,tau,work,lwork,info)
       CALL DSTERF(nn,diag,off_diag,info)

       lrii%cond_num=MAXVAL(ABS(diag))/MINVAL(ABS(diag))
       
       DEALLOCATE(diag,off_diag,smat,tau,work)
       CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    END DO 

    CALL neighbor_list_iterator_release(nl_iterator)


    CALL timestop(handle)

  END SUBROUTINE get_condition_number_of_overlap

! *****************************************************************************
!> \brief print recent information on optimization
!> \param opt_state state of the optimizer
!> \param lri_opt optimization environment
!> \param iunit ...
! *****************************************************************************
  SUBROUTINE print_optimization_update(opt_state,lri_opt,iunit)

    TYPE(opt_state_type)                     :: opt_state
    TYPE(lri_opt_type), POINTER              :: lri_opt
    INTEGER, INTENT(IN)                      :: iunit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_optimization_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n10

    n10 = MAX(opt_state%maxfun/100,1)

    IF ( opt_state%nf == 2 .AND. opt_state%state ==2 .AND. iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Initial value of function",T61,F20.10)') opt_state%f
    END IF
    IF ( MOD(opt_state%nf,n10) == 0 .AND. opt_state%nf > 1 .AND. iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
            INT(REAL(opt_state%nf,dp)/REAL(opt_state%maxfun,dp)*100._dp), opt_state%fopt
    ENDIF
    IF(lri_opt%use_condition_number) THEN
      IF ( MOD(opt_state%nf,n10) == 0 .AND. opt_state%nf > 1 .AND. iunit > 0 ) THEN
         WRITE(iunit,'(" POWELL| Recent value of function without condition nr.",T61,F20.10)') &
               lri_opt%rho_diff
      ENDIF
    ENDIF
  
  END SUBROUTINE  print_optimization_update

! *****************************************************************************
!> \brief write optimized LRI basis set to file
!> \param lri_env ...
!> \param dft_section ...
!> \param nkind ...
!> \param lri_opt ...
!> \param atomic_kind_set ...
! *****************************************************************************
  SUBROUTINE write_optimized_lri_basis(lri_env,dft_section,nkind,lri_opt,&
                                       atomic_kind_set)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(section_vals_type), POINTER         :: dft_section
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(lri_opt_type), POINTER              :: lri_opt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_optimized_lri_basis', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: cc_l, ikind, ipgf, iset, &
                                                ishell, nset, output_file
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc_orig
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: fbas
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY(fbas,gcc_orig,l,lmax,lmin,logger,npgf,nshell,print_key,zet)

    !*** do the printing 
    print_key => section_vals_get_subs_vals(dft_section,&
                                            "PRINT%OPTIMIZE_LRI_BASIS")
    logger => cp_get_default_logger()
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              dft_section,"PRINT%OPTIMIZE_LRI_BASIS"),&
              cp_p_file)) THEN
      output_file=cp_print_key_unit_nr(logger,dft_section,&
                                       "PRINT%OPTIMIZE_LRI_BASIS",&
                                       extension=".opt",&
                                       file_status="REPLACE",&
                                       file_action="WRITE",&
                                       file_form="FORMATTED") 

      IF(output_file>0) THEN

        filename = cp_print_key_generate_filename(logger,&
                   print_key,extension=".opt", &
                   my_local=.TRUE.) 
        
        DO ikind =1, nkind
          fbas  => lri_env%ri_basis(ikind)%gto_basis_set
          gcc_orig => lri_opt%ri_gcc_orig(ikind)%gcc_orig
          CALL get_gto_basis_set(gto_basis_set=fbas,&
                                 l=l, lmax=lmax, lmin=lmin,&
                                 npgf=npgf,nshell=nshell,&
                                 nset=nset,zet=zet)
          WRITE(output_file,'(T1,A2,T5,A)') TRIM(atomic_kind_set(ikind)%name),&
                                            TRIM(fbas%name)
          WRITE(output_file,'(T1,I4)') nset
          DO iset = 1, nset 
            WRITE(output_file,'(4(1X,I0))',advance='no') 2 ,lmin(iset), &
                                                         lmax(iset), npgf(iset)
            cc_l = 1
            DO ishell = 1,nshell(iset)
              IF(ishell /= nshell(iset)) THEN
                IF(l(ishell,iset)==l(ishell+1,iset)) THEN
                  cc_l = cc_l + 1
                ELSE
                  WRITE(output_file,'(1X,I0)',advance='no') cc_l
                  cc_l = 1
                ENDIF
              ELSE
                WRITE(output_file,'(1X,I0)') cc_l
              ENDIF
            ENDDO 
            DO ipgf = 1, npgf(iset)
              WRITE(output_file,'(F18.12)',advance='no') zet(ipgf,iset)
              DO ishell=1,nshell(iset)
                IF(ishell ==  nshell(iset)) THEN
                  WRITE(output_file,'(T5,F18.12)') gcc_orig(ipgf,ishell,iset)
                ELSE
                  WRITE(output_file,'(T5,F18.12)',advance='no') gcc_orig(ipgf,ishell,iset)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF
 
      CALL cp_print_key_finished_output(output_file,logger,dft_section,&
                        "PRINT%OPTIMIZE_LRI_BASIS")
    ENDIF

  END SUBROUTINE write_optimized_lri_basis

END MODULE lri_optimize_ri_basis
