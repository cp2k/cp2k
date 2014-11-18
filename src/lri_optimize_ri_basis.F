!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE lri_environment_methods,         ONLY: calculate_avec,&
                                             calculate_lri_integrals,&
                                             calculate_lri_overlap_aabb
  USE lri_environment_types,           ONLY: lri_basis_init,&
                                             lri_density_type,&
                                             lri_environment_type,&
                                             lri_int_rho_type,&
                                             lri_int_type,&
                                             lri_list_type,&
                                             lri_rhoab_type
  USE mathconstants,                   ONLY: pi
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE optimize_lri_basis(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_lri_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iunit, n10, nkind, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(opt_state_type)                     :: opt_state
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                lri_optbas_section

    NULLIFY(atomic_kind_set, dft_section, lri_density, lri_env, &
         lri_optbas_section, rho_struct)
    NULLIFY(input, logger, para_env, x)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,input=input,&
                    lri_env=lri_env,lri_density=lri_density,nkind=nkind,&
                    para_env=para_env,rho=rho_struct,error=error)

    ! density matrix
    CALL qs_rho_get(rho_struct, rho_ao=pmatrix, error=error)

    logger => cp_error_get_logger(error)
    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    lri_optbas_section => section_vals_get_subs_vals(input,&
                          "DFT%QS%OPTIMIZE_LRI_BASIS",error=error)
    iunit=cp_print_key_unit_nr(logger,input,"PRINT%PROGRAM_RUN_INFO",&
                               extension=".opt",error=error)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")')
    ENDIF

    ! *** initialization
    CALL init_optimization(lri_env,lri_optbas_section,&
                           opt_state,x,nkind,iunit,error)
    CALL calculate_lri_overlap_aabb(lri_env,qs_env,error)

    n10 = MAX(opt_state%maxfun/100,1)

    ! *** start optimize 
    opt_state%state = 0
    DO
      IF ( opt_state%state == 2 ) THEN
        CALL calc_lri_integrals_get_objective(lri_env,lri_density,qs_env,&
                                              opt_state,pmatrix,para_env,&
                                              atomic_kind_set,nkind,x,error)
      ENDIF

      IF ( opt_state%state == -1 ) EXIT

      ! *** ensure that exponents will be positive 
      x = SQRT(x)
      CALL powell_optimize (opt_state%nvar, x, opt_state)
      x= x**2._dp

      IF ( opt_state%nf == 2 .AND. opt_state%state ==2 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(/," POWELL| Initial value of function",T61,F20.10)') opt_state%f
      END IF
      IF ( MOD(opt_state%nf,n10) == 0 .AND. opt_state%nf > 1 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
              INT(REAL(opt_state%nf,dp)/REAL(opt_state%maxfun,dp)*100._dp), opt_state%fopt
      END IF
    ENDDO

    opt_state%state = 8

    x = SQRT(x)
    CALL powell_optimize (opt_state%nvar, x, opt_state)
    x= x**2._dp

    CALL write_optimized_lri_basis(lri_env,dft_section,nkind,x,&
         atomic_kind_set,error)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') opt_state%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,F20.10)') opt_state%fopt
      WRITE(iunit,'(/," Printed optimized lri basis set to file")')
    ENDIF

    CALL cp_print_key_finished_output(iunit,logger,input,&
                    "PRINT%PROGRAM_RUN_INFO", error=error)

    DEALLOCATE (x,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE optimize_lri_basis

! *****************************************************************************
!> \brief calculate objective only...
!> \param qs_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note  not called; testing only, to be removed later
! *****************************************************************************
  SUBROUTINE calc_objective_only(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_objective_only', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iunit, nkind
    REAL(KIND=dp)                            :: objf
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: input

    NULLIFY(atomic_kind_set, lri_density, lri_env, para_env, rho_struct)
    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,input=input,&
                    lri_env=lri_env,lri_density=lri_density,nkind=nkind,&
                    rho=rho_struct,para_env=para_env,error=error)
    CALL qs_rho_get(rho_struct, rho_ao=pmatrix, error=error)

    objf=0._dp 
    logger => cp_error_get_logger(error)
    iunit=cp_print_key_unit_nr(logger,input,"PRINT%PROGRAM_RUN_INFO",&
                               extension=".opt",error=error)
    CALL calculate_lri_overlap_aabb(lri_env,qs_env,error)
    CALL calculate_lri_integrals(lri_env,qs_env,calculate_forces=.FALSE.,error=error)
    CALL calculate_avec(lri_env,lri_density,qs_env,pmatrix,error=error)
    CALL calculate_objective(lri_env,lri_density,qs_env,pmatrix,para_env,objf,error)
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'("OBJF",T71,F20.15)') objf
    ENDIF
    CALL cp_print_key_finished_output(iunit,logger,input,&
                    "PRINT%PROGRAM_RUN_INFO", error=error)
  END SUBROUTINE calc_objective_only

! *****************************************************************************
!> \brief calculates the lri integrals and coefficients with the new exponents
!>        of the lri basis sets and calculates the objective function
!> \param lri_env ...
!> \param lri_optbas_section ...
!> \param opt_state ...
!> \param x parameters to be optimized, i.e. exponents of the lri basis set 
!> \param nkind ...
!> \param iunit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE init_optimization(lri_env,lri_optbas_section,opt_state,x,nkind,&
                               iunit,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(section_vals_type), POINTER         :: lri_optbas_section
    TYPE(opt_state_type)                     :: opt_state
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    INTEGER, INTENT(IN)                      :: nkind, iunit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_optimization', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, iset, n, nset, stat
    INTEGER, DIMENSION(:), POINTER           :: npgf
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    TYPE(gto_basis_set_type), POINTER        :: fbas

    NULLIFY(fbas,npgf,zet)
    failure = .FALSE.

    ALLOCATE(lri_env%ri_gcc_orig(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL section_vals_val_get(lri_optbas_section,"ACCURACY",&
                              r_val=opt_state%rhoend, error=error)
    CALL section_vals_val_get(lri_optbas_section,"STEP_SIZE",&
                              r_val=opt_state%rhobeg, error=error)
    CALL section_vals_val_get(lri_optbas_section,"MAX_FUN",&
                              i_val=opt_state%maxfun, error=error)

    opt_state%nvar =0
    opt_state%nf = 0
    opt_state%iprint = 1
    opt_state%unit  = iunit

    ! *** init exponents
    n = 0
    DO ikind=1,nkind
      fbas  => lri_env%ri_basis(ikind)%gto_basis_set
      CALL get_gto_basis_set(gto_basis_set=fbas,&
           npgf=npgf,nset=nset,zet=zet)
      DO iset =1,nset
        opt_state%nvar = opt_state%nvar + npgf(iset)
        CALL reallocate(x,1,opt_state%nvar)
        x(n+1:n+npgf(iset)) = zet(:,iset)
        n = n + npgf(iset)
      ENDDO
    ENDDO
    
    ! *** get the original gcc without normalization factor
    DO ikind=1,nkind
      fbas => lri_env%ri_basis(ikind)%gto_basis_set
      CALL get_original_gcc(lri_env%ri_gcc_orig(ikind)%gcc_orig,fbas,error)
    ENDDO

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
!> \brief calculates the lri integrals and coefficients with the new exponents
!>        of the lri basis sets and calculates the objective function
!> \param lri_env ...
!> \param lri_density ...
!> \param qs_env ...
!> \param opt_state ...
!> \param pmatrix density matrix
!> \param para_env ...
!> \param atomic_kind_set ...
!> \param nkind ...
!> \param x parameters to be optimized, i.e. exponents of the lri basis set 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calc_lri_integrals_get_objective(lri_env,lri_density,qs_env,&
                                              opt_state,pmatrix,para_env,&
                                              atomic_kind_set,nkind,x,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(opt_state_type)                     :: opt_state
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calc_lri_integrals_get_objective', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, ipgf, iset, n, nset
    INTEGER, DIMENSION(:), POINTER           :: npgf
    TYPE(gto_basis_set_type), POINTER        :: fbas

    NULLIFY(fbas,npgf)

    !*** sort the exponents and build new transformation matrices sphi
    n = 0
    DO ikind =1,nkind
      fbas  => lri_env%ri_basis(ikind)%gto_basis_set
      CALL get_gto_basis_set(gto_basis_set=fbas,npgf=npgf,nset=nset)
      DO iset =1,nset
        DO ipgf=1,npgf(iset)
          fbas%zet(ipgf,iset) = x(n+ipgf)
        ENDDO
        n = n + npgf(iset)
      ENDDO
      !build new sphi
      fbas%gcc = lri_env%ri_gcc_orig(ikind)%gcc_orig
      CALL init_orb_basis_set(fbas,error)
    ENDDO
    CALL lri_basis_init(lri_env,atomic_kind_set,error)
    CALL calculate_lri_integrals(lri_env,qs_env,calculate_forces=.FALSE.,error=error)
    CALL calculate_avec(lri_env,lri_density,qs_env,pmatrix,error=error)
    CALL calculate_objective(lri_env,lri_density,qs_env,pmatrix,para_env,&
                             opt_state%f,error)


  END SUBROUTINE calc_lri_integrals_get_objective

! *****************************************************************************
!> \brief calculates the objective function defined as integral of the square
!>        of rhoexact - rhofit, i.e. integral[(rhoexact-rhofit)**2]
!>        rhoexact is the exact pair density and rhofit the lri pair density
!> \param lri_env ...
!> \param lri_density ...
!> \param qs_env ...
!> \param pmatrix density matrix
!> \param para_env ...
!> \param fobj objective function
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calculate_objective(lri_env,lri_density,qs_env,pmatrix,para_env,&
                                 fobj,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(OUT)               :: fobj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
 
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
                                      *lri_env%bas_ovlp(ikind)%ri_ovlp(isgfa,jsgfa)
            ENDDO
          ENDDO
          IF(iatom /= jatom) THEN
            DO ksgfb=1,nfb
              DO lsgfb=1,nfb
                rhofit_sq = rhofit_sq + lrho%avec(nfa+ksgfb)*lrho%avec(nfa+lsgfb)&
                                       *lri_env%bas_ovlp(jkind)%ri_ovlp(ksgfb,lsgfb)
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

          fobj = fobj + obj_ab

        ENDDO

        CALL neighbor_list_iterator_release(nl_iterator)

      ENDDO
      CALL mp_sum(fobj,para_env%group) 

    ENDIF 

    CALL timestop(handle)

  END SUBROUTINE calculate_objective

! *****************************************************************************
!> \brief primitive Cartesian Gaussian functions are normalized. The normalization
!>        factor is included in the Gaussian contraction coefficients.
!>        Division by this factor to get the original gcc. 
!> \param gcc_orig ...
!> \param gto_basis_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE get_original_gcc(gcc_orig,gto_basis_set,error)
 
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc_orig
    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_original_gcc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, ishell, l, &
                                                maxpgf, maxshell, nset, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: expzet, gcca, prefac, zeta

    failure  = .FALSE.
    maxpgf   = SIZE(gto_basis_set%gcc,1)
    maxshell = SIZE(gto_basis_set%gcc,2)
    nset     = SIZE(gto_basis_set%gcc,3)

    ALLOCATE(gcc_orig(maxpgf,maxshell,nset),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
    DO iset=1,gto_basis_set%nset
       DO ishell=1,gto_basis_set%nshell(iset)
          l = gto_basis_set%l(ishell,iset)
          expzet = 0.25_dp*REAL(2*l + 3,dp)
          prefac = 2.0_dp**l*(2.0_dp/pi)**0.75_dp
          DO ipgf=1,gto_basis_set%npgf(iset)
             gcca = gto_basis_set%gcc(ipgf,ishell,iset)
             zeta = gto_basis_set%zet(ipgf,iset)
             gcc_orig(ipgf,ishell,iset) = gcca/(prefac*zeta**expzet)
          END DO
       END DO
    END DO

  END SUBROUTINE get_original_gcc

! *****************************************************************************
!> \brief write optimized LRI basis set to file
!> \param lri_env ...
!> \param dft_section ...
!> \param nkind ...
!> \param xopt ...
!> \param atomic_kind_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_optimized_lri_basis(lri_env,dft_section,nkind,xopt,&
                                       atomic_kind_set,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(section_vals_type), POINTER         :: dft_section
    INTEGER, INTENT(IN)                      :: nkind
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xopt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_optimized_lri_basis', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: cc_l, ikind, ipgf, iset, &
                                                ishell, n, nset, output_file
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc_orig
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: fbas
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY(fbas,gcc_orig,l,lmax,lmin,logger,npgf,nshell,print_key,zet)

    !*** sort the exponents
    n = 0
    DO ikind =1,nkind
      fbas  => lri_env%ri_basis(ikind)%gto_basis_set
      CALL get_gto_basis_set(gto_basis_set=fbas,&
           npgf=npgf,nset=nset,zet=zet)
      DO iset =1,nset
        DO ipgf=1,npgf(iset)
          zet(ipgf,iset) = xopt(n+ipgf)
        ENDDO
        n = n + npgf(iset)
      ENDDO
    ENDDO
   
    !*** do the printing 
    print_key => section_vals_get_subs_vals(dft_section,&
                                            "PRINT%OPTIMIZE_LRI_BASIS",&
                                            error=error)
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              dft_section,"PRINT%OPTIMIZE_LRI_BASIS",error=error),&
              cp_p_file)) THEN
      output_file=cp_print_key_unit_nr(logger,dft_section,&
                                       "PRINT%OPTIMIZE_LRI_BASIS",&
                                       extension=".opt",&
                                       file_status="REPLACE",&
                                       file_action="WRITE",&
                                       file_form="FORMATTED",& 
                                       error=error)

      IF(output_file>0) THEN

        filename = cp_print_key_generate_filename(logger,&
                   print_key,extension=".opt", &
                   my_local=.TRUE.,error=error) 
        
        DO ikind =1, nkind
          fbas  => lri_env%ri_basis(ikind)%gto_basis_set
          gcc_orig => lri_env%ri_gcc_orig(ikind)%gcc_orig
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
                        "PRINT%OPTIMIZE_LRI_BASIS", error=error)
    ENDIF

  END SUBROUTINE write_optimized_lri_basis

END MODULE lri_optimize_ri_basis
