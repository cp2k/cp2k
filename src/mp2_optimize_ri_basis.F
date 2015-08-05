!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to optimize the RI-MP2 basis. Only exponents of  
!>        non-contracted auxiliary basis basis are optimized. The derivative
!>        of the MP2 energy with respect to the exponents of the basis
!>        are calculated numerically.
!> \par History
!>      08.2013 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
MODULE mp2_optimize_ri_basis
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE basis_set_container_types,       ONLY: add_basis_set_to_container,&
                                             remove_basis_from_container
  USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
! 
  USE hfx_libint_wrapper_types,        ONLY: build_eri_size
  USE hfx_types,                       ONLY: hfx_basis_info_type,&
                                             hfx_basis_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: default_output_unit
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_comm_split_direct,&
                                             mp_sum
  USE mp2_direct_method,               ONLY: mp2_canonical_direct_single_batch
  USE mp2_ri_libint,                   ONLY: libint_ri_mp2,&
                                             read_RI_basis_set,&
                                             release_RI_basis_set
  USE mp2_types,                       ONLY: mp2_biel_type,&
                                             mp2_type
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset,&
                                             nso
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_optimize_ri_basis'

  PUBLIC :: optimize_ri_basis_main

  CONTAINS

! *****************************************************************************
!> \brief optimize RI-MP2 basis set
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_ex ...
!> \param Emp2_S ...
!> \param Emp2_T ...
!> \param dimen ...
!> \param natom ...
!> \param homo ...
!> \param mp2_biel ...
!> \param mp2_env ...
!> \param C ...
!> \param Auto ...
!> \param kind_of ...
!> \param qs_env ...
!> \param rho ...
!> \param hfx_sections ...
!> \param para_env ...
!> \param unit_nr ...
!> \param error ...
!> \param homo_beta ...
!> \param C_beta ...
!> \param Auto_beta ...
!> \author Mauro Del Ben
! *****************************************************************************
 SUBROUTINE optimize_ri_basis_main(Emp2,Emp2_Cou,Emp2_ex,Emp2_S,Emp2_T,dimen,natom,homo, &
                                     mp2_biel,mp2_env,C,Auto,kind_of,&
                                     qs_env,rho,hfx_sections,para_env, &
                                     unit_nr,error,homo_beta,C_beta,Auto_beta)

    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_ex, &
                                                Emp2_S, Emp2_T
    INTEGER                                  :: dimen, natom, homo
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, OPTIONAL                        :: homo_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :), OPTIONAL              :: C_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), OPTIONAL                 :: Auto_beta

    CHARACTER(len=*), PARAMETER :: routineN = 'optimize_ri_basis_main', &
      routineP = moduleN//':'//routineN

    INTEGER :: color_sub, comm_sub, dimen_RI, elements_ij_proc, handle, i, &
      iiter, ikind, ipgf, iset, ishell, j, local_unit_nr, max_l_am, &
      max_num_iter, ndof, nkind, number_groups, stat, virtual, virtual_beta
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_list_proc, index_table_RI
    LOGICAL                                  :: failure, hess_up, &
                                                open_shell_case, &
                                                reset_boundary
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: basis_was_assoc
    REAL(KIND=dp) :: DI, DI_new, DRI, DRI_new, Emp2_AA, Emp2_AA_Cou, &
      Emp2_AA_ex, Emp2_AB, Emp2_AB_Cou, Emp2_AB_ex, Emp2_BB, Emp2_BB_Cou, &
      Emp2_BB_ex, Emp2_RI, Emp2_RI_new, eps_DI_rel, eps_DRI, eps_step, expon, &
      fac, fad, fae, reset_edge, sumdg, sumxi
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: deriv, dg, g, hdg, lower_B, &
                                                max_dev, max_rel_dev, p, &
                                                pnew, xi
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: exp_limits, hessin
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: Integ_MP2, Integ_MP2_AA, &
                                                Integ_MP2_AB, Integ_MP2_BB
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type)                      :: error_sub
    TYPE(cp_logger_type), POINTER            :: logger, logger_sub
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    TYPE(gto_basis_set_type), POINTER        :: ri_aux_basis
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0, RI_basis_parameter
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)

    open_shell_case=.FALSE.
    IF(PRESENT(homo_beta).AND.PRESENT(C_beta).AND.PRESENT(Auto_beta)) open_shell_case=.TRUE.

    virtual=dimen-homo

    eps_DRI=mp2_env%ri_opt_param%DRI
    eps_DI_rel=mp2_env%ri_opt_param%DI_rel
    eps_step=mp2_env%ri_opt_param%eps_step
    max_num_iter=mp2_env%ri_opt_param%max_num_iter

    ! calculate the ERI's over molecular integrals
    Emp2=0.0_dp
    Emp2_Cou=0.0_dp
    Emp2_ex=0.0_dp
    Emp2_S=0.0_dp
    Emp2_T=0.0_dp
    IF(open_shell_case) THEN
      ! open shell case
      virtual_beta=dimen-homo_beta

      ! alpha-aplha case
      Emp2_AA=0.0_dp
      Emp2_AA_Cou=0.0_dp
      Emp2_AA_ex=0.0_dp
      CALL calc_elem_ij_proc(homo,homo,para_env,elements_ij_proc,ij_list_proc,error)
      CALL mp2_canonical_direct_single_batch(Emp2_AA,Emp2_AA_Cou,Emp2_AA_ex,mp2_env,qs_env,rho,hfx_sections,para_env,&
                                             mp2_biel,dimen,C,Auto,0,homo,homo,&
                                             elements_ij_proc,ij_list_proc,homo,0,&
                                             Integ_MP2=Integ_MP2_AA,error=error)
      CALL mp_sum(Emp2_AA_Cou,para_env%group)
      CALL mp_sum(Emp2_AA_Ex,para_env%group)
      CALL mp_sum(Emp2_AA,para_env%group)
      DEALLOCATE(ij_list_proc)

      ! beta-beta case
      Emp2_BB=0.0_dp
      Emp2_BB_Cou=0.0_dp
      Emp2_BB_ex=0.0_dp
      CALL calc_elem_ij_proc(homo_beta,homo_beta,para_env,elements_ij_proc,ij_list_proc,error)
      CALL mp2_canonical_direct_single_batch(Emp2_BB,Emp2_BB_Cou,Emp2_BB_ex,mp2_env,qs_env,rho,hfx_sections,para_env,&
                                             mp2_biel,dimen,C_beta,Auto_beta,0,homo_beta,homo_beta,&
                                             elements_ij_proc,ij_list_proc,homo_beta,0,&
                                             Integ_MP2=Integ_MP2_BB,error=error)
      CALL mp_sum(Emp2_BB_Cou,para_env%group)
      CALL mp_sum(Emp2_BB_Ex,para_env%group)
      CALL mp_sum(Emp2_BB,para_env%group)
      DEALLOCATE(ij_list_proc)

      ! aplha-beta case
      Emp2_AB=0.0_dp
      Emp2_AB_Cou=0.0_dp
      Emp2_AB_ex=0.0_dp
      CALL calc_elem_ij_proc(homo,homo_beta,para_env,elements_ij_proc,ij_list_proc,error)
      CALL mp2_canonical_direct_single_batch(Emp2_AB,Emp2_AB_Cou,Emp2_AB_ex,mp2_env,qs_env,rho,hfx_sections,para_env,&
                                             mp2_biel,dimen,C,Auto,0,homo,homo,&
                                             elements_ij_proc, ij_list_proc,homo_beta,0,&
                                             homo_beta,C_beta,Auto_beta,Integ_MP2=Integ_MP2_AB,error=error)
      CALL mp_sum(Emp2_AB_Cou,para_env%group)
      CALL mp_sum(Emp2_AB_Ex,para_env%group)
      CALL mp_sum(Emp2_AB,para_env%group)
      DEALLOCATE(ij_list_proc)

      ! IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Alpha = ', Emp2_AA
      ! IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'MP2 Energy Beta-Beta   = ', Emp2_BB
      ! IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Beta  = ', Emp2_AB

      Emp2=Emp2_AA+Emp2_BB+Emp2_AB*2.0_dp !+Emp2_BA
      Emp2_Cou=Emp2_AA_Cou+Emp2_BB_Cou+Emp2_AB_Cou*2.0_dp !+Emp2_BA
      Emp2_ex=Emp2_AA_ex+Emp2_BB_ex+Emp2_AB_ex*2.0_dp !+Emp2_BA

      Emp2_S=Emp2_AB*2.0_dp
      Emp2_T=Emp2_AA+Emp2_BB

      ! Replicate the MO-ERI's over all processes
      CALL mp_sum(Integ_MP2_AA,para_env%group)
      CALL mp_sum(Integ_MP2_BB,para_env%group)
      CALL mp_sum(Integ_MP2_AB,para_env%group)

    ELSE
      ! close shell case
      CALL calc_elem_ij_proc(homo,homo,para_env,elements_ij_proc,ij_list_proc,error)
      CALL mp2_canonical_direct_single_batch(Emp2,Emp2_Cou,Emp2_ex,mp2_env,qs_env,rho,hfx_sections,para_env,&
                                             mp2_biel,dimen,C,Auto,0,homo,homo,&
                                             elements_ij_proc,ij_list_proc,homo,0,&
                                             Integ_MP2=Integ_MP2,error=error)
      CALL mp_sum(Emp2_Cou,para_env%group)
      CALL mp_sum(Emp2_Ex,para_env%group)
      CALL mp_sum(Emp2,para_env%group)
      DEALLOCATE(ij_list_proc)

      ! Replicate the MO-ERI's over all processes
      CALL mp_sum(Integ_MP2,para_env%group)

    END IF

    ! create the para_env_sub
    number_groups=para_env%num_pe/mp2_env%mp2_num_proc
    color_sub=para_env%mepos/mp2_env%mp2_num_proc
    CALL mp_comm_split_direct(para_env%group,comm_sub,color_sub)
    NULLIFY(para_env_sub)
    CALL cp_para_env_create(para_env_sub,comm_sub,error=error)

    IF (para_env%mepos==para_env%source) THEN
       local_unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
    ELSE
       local_unit_nr=default_output_unit
    ENDIF
    NULLIFY(logger_sub)
    CALL cp_logger_create(logger_sub,para_env=para_env_sub,&
               default_global_unit_nr=local_unit_nr, close_global_unit_on_dealloc=.FALSE.)
    CALL cp_logger_set(logger_sub,local_filename="opt_RI_basis_localLog")
    CALL cp_error_init(error_sub, logger=logger_sub)

    CALL generate_RI_init_basis(qs_env,mp2_env,nkind,max_rel_dev,basis_was_assoc,error)

    CALL read_RI_basis_set(qs_env,RI_basis_parameter,RI_basis_info,&
                           natom,nkind,kind_of,index_table_RI,dimen_RI,&
                           basis_S0,error) 

    ndof=0
    max_l_am=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ndof=ndof+1
        max_l_am=MAX(max_l_am,MAXVAL(RI_basis_parameter(ikind)%lmax))
      END DO
    END DO
 
    ALLOCATE(exp_limits(2,nkind))
    exp_limits(1,:)= HUGE(0)
    exp_limits(2,:)=-HUGE(0) 
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        expon=RI_basis_parameter(ikind)%zet(1,iset)
        IF(expon<=exp_limits(1,ikind)) exp_limits(1,ikind)=expon
        IF(expon>=exp_limits(2,ikind)) exp_limits(2,ikind)=expon
      END DO
      IF(basis_was_assoc(ikind)) THEN
        exp_limits(1,ikind)=exp_limits(1,ikind)*0.5_dp
        exp_limits(2,ikind)=exp_limits(2,ikind)*1.5_dp
      ELSE
        exp_limits(1,ikind)=exp_limits(1,ikind)*0.8_dp*0.5_dp
        exp_limits(2,ikind)=exp_limits(2,ikind)*1.2_dp*1.5_dp
      END IF
    END DO
    DEALLOCATE(basis_was_assoc)

    ! check if the max angular momentum exceed the libint one
    IF(max_l_am>build_eri_size) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The angular momentum needed exceeds the value assumed when configuring libint.")
    END IF

    ! Allocate stuff
    ALLOCATE(p(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    p=0.0_dp
    ALLOCATE(xi(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    xi=0.0_dp
    ALLOCATE(g(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    g=0.0_dp
    ALLOCATE(dg(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dg=0.0_dp
    ALLOCATE(hdg(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    hdg=0.0_dp
    ALLOCATE(pnew(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    pnew=0.0_dp
    ALLOCATE(deriv(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    deriv=0.0_dp
    ALLOCATE(hessin(ndof,ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    hessin=0.0_dp
    DO i=1, ndof 
      hessin(i,i)=1.0_dp
    END DO

    ! initialize trasformation function
    ALLOCATE(lower_B(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    lower_B=0.0_dp
    ALLOCATE(max_dev(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    max_dev=0.0_dp

    ! Initialize the transformation function
    CALL init_transf(nkind,RI_basis_parameter,lower_B,max_dev,max_rel_dev)

    ! get the atomic kind set for writing the basis
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,error=error)

    ! Calculate RI-MO-ERI's
    CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                          Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                          qs_env,natom,dimen,dimen_RI,homo,virtual,&
                          kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                          RI_basis_parameter,RI_basis_info,basis_S0,&
                          open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                          .TRUE.,error)

    ! ! Calculate function (DI) derivatives with respect to the RI basis exponent
    CALL calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps_step,&
                              qs_env,nkind,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                              RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                              para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                              p,lower_B,max_dev,&
                              deriv,error,error_sub)

    g(:)=deriv
    xi(:)=-g

    reset_edge=1.5_dp
    DO iiter=1, max_num_iter
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,I5)') 'OPTIMIZATION STEP NUMBER',iiter

      ! perform step
      pnew(:)=p+xi
      CALL p2basis(nkind,RI_basis_parameter,lower_B,max_dev,pnew)

      ! check if we have to reset boundaries 
      reset_boundary=.FALSE.
      ! nreset=0
      i=0
      DO ikind=1, nkind
        DO iset=1, RI_basis_parameter(ikind)%nset
          i=i+1
          CALL transf_val(lower_B(i), max_dev(i), pnew(i), expon)
          IF(ABS(pnew(i))>reset_edge.OR.expon<exp_limits(1,ikind).OR.expon>exp_limits(2,ikind)) THEN
            reset_boundary=.TRUE.
            EXIT
          END IF
        END DO
      END DO
      ! IF(nreset>1) reset_boundary=.TRUE.

      IF(reset_boundary) THEN
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'RESET BASIS: one of the exponent hits the boundary'
        CALL reset_basis(nkind,ndof,RI_basis_parameter,reset_edge,&
                         pnew,lower_B,max_dev,max_rel_dev,exp_limits,error)
        p(:)=pnew
        xi=0.0_dp
        g=0.0_dp
        dg=0.0_dp
        hdg=0.0_dp
        pnew=0.0_dp
        hessin=0.0_dp
        DO i=1, ndof
          hessin(i,i)=1.0_dp
        END DO
        deriv=0.0_dp
        ! Calculate RI-MO-ERI's
        CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                              qs_env,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                              RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                              .FALSE.,error)
        ! ! Calculate function (DI) derivatives with respect to the RI basis exponent
        CALL calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,DI,&
                                  Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps_step,&
                                  qs_env,nkind,natom,dimen,dimen_RI,homo,virtual,&
                                  kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                                  RI_basis_parameter,RI_basis_info,basis_S0,&
                                  open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                                  para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                                  p,lower_B,max_dev,&
                                  deriv,error,error_sub)

        g(:)=deriv
        xi(:)=-g
        pnew(:)=p+xi
        CALL p2basis(nkind,RI_basis_parameter,lower_B,max_dev,pnew)
      END IF

      ! calculate energy at the new point
      CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI_new,DRI_new,DI_new,&
                            Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                            qs_env,natom,dimen,dimen_RI,homo,virtual,&
                            kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                            RI_basis_parameter,RI_basis_info,basis_S0,&
                            open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                            .FALSE.,error)
 
      ! update energy and direction
      DI=DI_new
      xi(:)=pnew-p
      p(:)=pnew

      ! check for convergence
      IF (unit_nr>0) THEN
        WRITE(unit_nr,*)
        DO ikind=1, nkind
          CALL get_qs_kind(qs_kind=qs_kind_set(ikind),basis_set=ri_aux_basis,&
               basis_type="RI_AUX",error=error)
          WRITE(unit_nr,'(T3,A,A)') atomic_kind_set(ikind)%element_symbol,'   RI_opt_basis'
          WRITE(unit_nr,'(T3,I3)') RI_basis_parameter(ikind)%nset
          DO iset=1, RI_basis_parameter(ikind)%nset
            WRITE(unit_nr,'(T3,10I4)') iset,&
                                       RI_basis_parameter(ikind)%lmin(iset),&
                                       RI_basis_parameter(ikind)%lmax(iset),&
                                       RI_basis_parameter(ikind)%npgf(iset),&
                                       (1,ishell=1,RI_basis_parameter(ikind)%nshell(iset))
            DO ipgf=1,RI_basis_parameter(ikind)%npgf(iset)
              WRITE(unit_nr,'(T3,10F16.10)') RI_basis_parameter(ikind)%zet(ipgf,iset),&
                                             (ri_aux_basis%gcc(ipgf,ishell,iset),&
                                             ishell=1,ri_aux_basis%nshell(iset))
            END DO                            
          END DO
          WRITE(unit_nr,*)
        END DO
        WRITE(unit_nr,*)
      END IF
      IF(DI/ABS(Emp2)<=eps_DI_rel.AND.ABS(DRI_new)<=eps_DRI) THEN
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'OPTIMIZATION CONVERGED'
        IF (unit_nr>0) WRITE(unit_nr,*)
        EXIT
      END IF

      ! calculate gradients
      CALL calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,DI,&
                                Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps_step,&
                                qs_env,nkind,natom,dimen,dimen_RI,homo,virtual,&
                                kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                                RI_basis_parameter,RI_basis_info,basis_S0,&
                                open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                                para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                                p,lower_B,max_dev,&
                                deriv,error,error_sub)

      ! g is the vector containing the old gradient
      dg(:)=deriv-g
      g(:)=deriv
      hdg(:)=MATMUL(hessin,dg)

      fac=SUM(dg*xi)
      fae=SUM(dg*hdg)
      sumdg=SUM(dg*dg)
      sumxi=SUM(xi*xi)
      
      hess_up=.TRUE.
      IF(fac**2>sumdg*sumxi*3.0E-8_dp) THEN
        fac=1.0_dp/fac
        fad=1.0_dp/fae
        dg(:)=fac*xi-fad*hdg
        DO i=1, ndof
          DO j=1, ndof
            hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)&
                                   -fad*hdg(i)*hdg(j)&
                                   +fae*dg(i)*dg(j)
          END DO
        END DO
      ELSE
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Skip Hessian Update'   
        hess_up=.FALSE.
      END IF

      !XXXXXXXXXXXXX
      ! DO i=1, ndof
      !   IF(g(i)<0.0D+00) THEN
      !     xi(i)=MAX(0.05,ABS(g(i)))
      !   ELSE
      !     xi(i)=-MAX(0.05,ABS(g(i)))
      !   END IF
      ! END DO
      !XXXXXXXXXXXXX

      ! new direction
      xi(:)=-MATMUL(hessin,g)

    END DO

    IF(.NOT.(DI/ABS(Emp2)<=eps_DI_rel.AND.ABS(DRI_new)<=eps_DRI)) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,I5,A)') 'OPTIMIZATION NOT CONVERGED IN',max_num_iter,' STEPS.'
      IF (unit_nr>0) WRITE(unit_nr,*)
    END IF

    DEALLOCATE(max_rel_dev)

    DEALLOCATE(p)
    DEALLOCATE(xi)
    DEALLOCATE(g)
    DEALLOCATE(pnew)
    DEALLOCATE(dg)  
    DEALLOCATE(hdg)  
    DEALLOCATE(deriv)
    DEALLOCATE(Hessin)
    DEALLOCATE(lower_B)
    DEALLOCATE(max_dev)
    DEALLOCATE(exp_limits)

    IF(open_shell_case) THEN
      DEALLOCATE(Integ_MP2_AA)
      DEALLOCATE(Integ_MP2_BB)
      DEALLOCATE(Integ_MP2_AB)
    ELSE
      DEALLOCATE(Integ_MP2)
    END IF
    DEALLOCATE(index_table_RI)

    ! Release RI basis set
    CALL release_RI_basis_set(RI_basis_parameter,basis_S0,error)
 
    CALL cp_para_env_release(para_env_sub,error)
    CALL cp_error_dealloc_ref(error_sub,error=error)
    CALL cp_logger_release(logger_sub)

    CALL timestop(handle)

  END SUBROUTINE optimize_ri_basis_main

! *****************************************************************************
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_AA ...
!> \param Emp2_BB ...
!> \param Emp2_AB ...
!> \param DI_ref ...
!> \param Integ_MP2 ...
!> \param Integ_MP2_AA ...
!> \param Integ_MP2_BB ...
!> \param Integ_MP2_AB ...
!> \param eps ...
!> \param qs_env ...
!> \param nkind ...
!> \param natom ...
!> \param dimen ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param kind_of ...
!> \param index_table_RI ...
!> \param mp2_biel ...
!> \param mp2_env ...
!> \param Auto ...
!> \param C ...
!> \param RI_basis_parameter ...
!> \param RI_basis_info ...
!> \param basis_S0 ...
!> \param open_shell_case ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param Auto_beta ...
!> \param C_beta ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param number_groups ...
!> \param color_sub ...
!> \param unit_nr ...
!> \param p ...
!> \param lower_B ...
!> \param max_dev ...
!> \param deriv ...
!> \param error ...
!> \param error_sub ...
! *****************************************************************************
  SUBROUTINE calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,DI_ref,&
                                  Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps,&
                                  qs_env,nkind,natom,dimen,dimen_RI,homo,virtual,&
                                  kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                                  RI_basis_parameter,RI_basis_info,basis_S0,&
                                  open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                                  para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                                  p,lower_B,max_dev,&
                                  deriv,error,error_sub)
    REAL(KIND=dp)                            :: Emp2, Emp2_AA, Emp2_BB, &
                                                Emp2_AB, DI_ref
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: Integ_MP2, Integ_MP2_AA, &
                                                Integ_MP2_BB, Integ_MP2_AB
    REAL(KIND=dp)                            :: eps
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER                                  :: nkind, natom, dimen, &
                                                dimen_RI, homo, virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: index_table_RI
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0
    LOGICAL                                  :: open_shell_case
    INTEGER                                  :: homo_beta, virtual_beta
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C_beta
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: number_groups, color_sub, &
                                                unit_nr
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p, lower_B, max_dev, deriv
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_energy_func_der', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ideriv, ikind, iset, &
                                                nseta
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: DI, DRI, Emp2_RI, &
                                                new_basis_val, orig_basis_val
    REAL(KIND=dp), VOLATILE                  :: step, temp

    CALL timeset(routineN,handle)
    failure=.FALSE.

    step=eps

    ! cycle over the RI basis set exponent
    deriv=0.0_dp
    ideriv=0
    DO ikind=1, nkind
      nseta=RI_basis_parameter(ikind)%nset
      DO iset=1, nseta
        ! for now only uncontracted aux basis set
        ideriv=ideriv+1
        IF (MOD(ideriv,number_groups)/=color_sub) CYCLE

        ! ! calculate the numerical derivative
        ! ! The eps is the relative change of the exponent for the
        ! ! calculation of the numerical derivative
        ! CPPostcondition(RI_basis_parameter(ikind)%npgf(iset)==1,cp_failure_level,routineP,error,failure)
        ! step=eps*RI_basis_parameter(ikind)%zet(1,iset)
        ! temp=RI_basis_parameter(ikind)%zet(1,iset)+step
        ! step=temp-RI_basis_parameter(ikind)%zet(1,iset)
        ! RI_basis_parameter(ikind)%zet(1,iset)=RI_basis_parameter(ikind)%zet(1,iset)+step

        ! in the new case eps is just the step length for calculating the numerical derivative

        CPPostcondition(RI_basis_parameter(ikind)%npgf(iset)==1,cp_failure_level,routineP,error,failure)
        orig_basis_val=RI_basis_parameter(ikind)%zet(1,iset)
        temp=p(ideriv)+step
        CALL transf_val(lower_B(ideriv), max_dev(ideriv), temp, new_basis_val)
        RI_basis_parameter(ikind)%zet(1,iset)=new_basis_val
     
        CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                              qs_env,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                              RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                              para_env_sub,unit_nr,.TRUE.,error_sub)

        RI_basis_parameter(ikind)%zet(1,iset)=orig_basis_val

        IF(para_env_sub%mepos==0) THEN
          temp=EXP(DI)
          temp=temp/EXP(DI_ref)
          deriv(ideriv)=LOG(temp)/step
        END IF

      END DO
    END DO

    CALL mp_sum(deriv,para_env%group)

    CALL timestop(handle)

  END SUBROUTINE 

! *****************************************************************************
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_AA ...
!> \param Emp2_BB ...
!> \param Emp2_AB ...
!> \param Emp2_RI ...
!> \param DRI ...
!> \param DI ...
!> \param Integ_MP2 ...
!> \param Integ_MP2_AA ...
!> \param Integ_MP2_BB ...
!> \param Integ_MP2_AB ...
!> \param qs_env ...
!> \param natom ...
!> \param dimen ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param kind_of ...
!> \param index_table_RI ...
!> \param mp2_biel ...
!> \param mp2_env ...
!> \param Auto ...
!> \param C ...
!> \param RI_basis_parameter ...
!> \param RI_basis_info ...
!> \param basis_S0 ...
!> \param open_shell_case ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param Auto_beta ...
!> \param C_beta ...
!> \param para_env ...
!> \param unit_nr ...
!> \param no_write ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                              qs_env,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,mp2_biel,mp2_env,Auto,C,&
                              RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                              no_write,error)
    REAL(KIND=dp)                            :: Emp2, Emp2_AA, Emp2_BB, &
                                                Emp2_AB, Emp2_RI, DRI, DI
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: Integ_MP2, Integ_MP2_AA, &
                                                Integ_MP2_BB, Integ_MP2_AB
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER                                  :: natom, dimen, dimen_RI, homo, &
                                                virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: index_table_RI
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0
    LOGICAL                                  :: open_shell_case
    INTEGER                                  :: homo_beta, virtual_beta
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C_beta
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: unit_nr
    LOGICAL                                  :: no_write
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_energy_func', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: DI_AA, DI_AB, DI_BB, DRI_AA, &
                                                DRI_AB, DRI_BB, Emp2_RI_AA, &
                                                Emp2_RI_AB, Emp2_RI_BB
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Lai, Lai_beta

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CALL libint_ri_mp2(dimen,dimen_RI,homo,natom,mp2_biel,mp2_env,C,&
                       kind_of,RI_basis_parameter,RI_basis_info,basis_S0,index_table_RI,&
                       qs_env,para_env,Lai,error)
    IF(open_shell_case) THEN
      CALL libint_ri_mp2(dimen,dimen_RI,homo_beta,natom,mp2_biel,mp2_env,C_beta,&
                         kind_of,RI_basis_parameter,RI_basis_info,basis_S0,index_table_RI,&
                         qs_env,para_env,Lai_beta,error)
    END IF

    ! Contract integrals into energy
    IF(open_shell_case) THEN
      ! alpha-alpha
      CALL contract_integrals(DI_AA,Emp2_RI_AA,DRI_AA,Emp2_AA,homo,homo,virtual,virtual,&
                              1.0_dp,0.5_dp,.TRUE.,&
                              Auto,Auto,Integ_MP2_AA,&
                              Lai,Lai,para_env,error)

      ! beta-beta
      CALL contract_integrals(DI_BB,Emp2_RI_BB,DRI_BB,Emp2_BB,homo_beta,homo_beta,virtual_beta,virtual_beta,&
                              1.0_dp,0.5_dp,.TRUE.,&
                              Auto_beta,Auto_beta,Integ_MP2_BB,&
                              Lai_beta,Lai_beta,para_env,error)

      ! alpha-beta
      CALL contract_integrals(DI_AB,Emp2_RI_AB,DRI_AB,Emp2_AB*2.0_dp,homo,homo_beta,virtual,virtual_beta,&
                              1.0_dp,1.0_dp,.FALSE.,&
                              Auto,Auto_beta,Integ_MP2_AB,&
                              Lai,Lai_beta,para_env,error)

      Emp2_RI=Emp2_RI_AA+Emp2_RI_BB+Emp2_RI_AB
      DRI=DRI_AA+DRI_BB+DRI_AB
      DI=DI_AA+DI_BB+DI_AB
    ELSE
      CALL contract_integrals(DI,Emp2_RI,DRI,Emp2,homo,homo,virtual,virtual,&
                              2.0_dp,1.0_dp,.TRUE.,&
                              Auto,Auto,Integ_MP2,&
                              Lai,Lai,para_env,error)
    END IF

    IF(.NOT.no_write) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)')
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'Emp2 =   ', Emp2
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'Emp2-RI =', Emp2_RI
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,ES25.10)') 'DRI =    ', DRI
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,ES25.10)') 'DI =     ', DI
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,ES25.10)') 'DI/|Emp2| =     ', DI/ABS(Emp2)
    END IF

    DEALLOCATE(Lai)
    IF(open_shell_case) DEALLOCATE(Lai_beta)
 
    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param nkind ...
!> \param RI_basis_parameter ...
!> \param lower_B ...
!> \param max_dev ...
!> \param max_rel_dev ...
! *****************************************************************************
  SUBROUTINE init_transf(nkind,RI_basis_parameter,lower_B,max_dev,max_rel_dev)
    INTEGER                                  :: nkind
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: lower_B, max_dev, max_rel_dev

    INTEGER                                  :: ikind, ipos, iset

    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        lower_B(ipos)=RI_basis_parameter(ikind)%zet(1,iset)*(1.0_dp-max_rel_dev(ipos))
        max_dev(ipos)=RI_basis_parameter(ikind)%zet(1,iset)*2.0_dp*max_rel_dev(ipos)
      END DO
    END DO

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param nkind ...
!> \param RI_basis_parameter ...
!> \param Lower_B ...
!> \param max_dev ...
!> \param p ...
! *****************************************************************************
  SUBROUTINE p2basis(nkind,RI_basis_parameter,Lower_B,max_dev,p)
    INTEGER                                  :: nkind
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Lower_B, max_dev, p

    INTEGER                                  :: ikind, ipos, iset
    REAL(KIND=dp)                            :: valout

    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        CALL transf_val(lower_B(ipos), max_dev(ipos), p(ipos), valout)
        RI_basis_parameter(ikind)%zet(1,iset)=valout
      END DO
    END DO

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param nkind ...
!> \param ndof ...
!> \param RI_basis_parameter ...
!> \param reset_edge ...
!> \param pnew ...
!> \param lower_B ...
!> \param max_dev ...
!> \param max_rel_dev ...
!> \param exp_limits ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE reset_basis(nkind,ndof,RI_basis_parameter,reset_edge,&
                         pnew,lower_B,max_dev,max_rel_dev,exp_limits,error)
    INTEGER                                  :: nkind, ndof
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    REAL(KIND=dp)                            :: reset_edge
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pnew, Lower_B, max_dev, &
                                                max_rel_dev
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: exp_limits
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reset_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: am_max, iexpo, ikind, ipos, &
                                                ipos_p, iset, la, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nf_per_l
    INTEGER, DIMENSION(ndof)                 :: change_expo
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: has_to_be_changed
    REAL(KIND=dp)                            :: expo, geom_fact, pmax
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: max_min_exp_per_l
    REAL(KIND=dp), DIMENSION(ndof)           :: new_expo, old_expo, &
                                                old_lower_B, old_max_dev, &
                                                old_max_rel_dev, old_pnew

    failure=.FALSE.

    ! make a copy of the original parameters
    old_pnew        = pnew
    old_lower_B     = lower_B
    old_max_dev     = max_dev
    old_max_rel_dev = max_rel_dev
    old_expo=0.0_dp
    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        old_expo(ipos)=RI_basis_parameter(ikind)%zet(1,iset)
      END DO
    END DO

    pnew = 0.0_dp
    lower_B = 0.0_dp
    max_dev = 0.0_dp
    max_rel_dev =0.0_dp

    change_expo=0

    new_expo=0.0_dp
    ipos=0
    ipos_p=0
    DO ikind=1, nkind
      am_max=MAXVAL(RI_basis_parameter(ikind)%lmax(:))
      ALLOCATE(nf_per_l(0:am_max),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      nf_per_l=0
      ALLOCATE(max_min_exp_per_l(2,0:am_max),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      max_min_exp_per_l(1,:)=HUGE(0)
      max_min_exp_per_l(2,:)=-HUGE(0)

      DO iset=1, RI_basis_parameter(ikind)%nset
        la=RI_basis_parameter(ikind)%lmax(iset)
        expo=RI_basis_parameter(ikind)%zet(1,iset)
        nf_per_l(la)=nf_per_l(la)+1
        IF(expo<=max_min_exp_per_l(1,la)) max_min_exp_per_l(1,la)=expo
        IF(expo>=max_min_exp_per_l(2,la)) max_min_exp_per_l(2,la)=expo
      END DO

      max_min_exp_per_l(1,la)=MAX(max_min_exp_per_l(1,la),exp_limits(1,ikind))
      max_min_exp_per_l(2,la)=MIN(max_min_exp_per_l(2,la),exp_limits(2,ikind))

      ! always s exponents as maximum and minimu
      ! expo=MAXVAL(max_min_exp_per_l(2,:))
      ! max_min_exp_per_l(2,0)=expo
      ! expo=MINVAL(max_min_exp_per_l(1,:))
      ! max_min_exp_per_l(1,0)=expo

      ALLOCATE(has_to_be_changed(0:am_max),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      has_to_be_changed=.FALSE.
      DO la=0, am_max
        pmax=-HUGE(0)
        DO iexpo=1, nf_per_l(la)
          ipos_p=ipos_p+1
          IF(ABS(old_pnew(ipos_p))>=pmax) pmax=ABS(old_pnew(ipos_p))
          ! check if any of the exponents go out of range
          CALL transf_val(old_lower_B(ipos_p),old_max_dev(ipos_p),old_pnew(ipos_p),expo)
          IF(expo<exp_limits(1,ikind).OR.expo>exp_limits(2,ikind)) has_to_be_changed(la)=.TRUE.
        END DO
        IF(pmax>reset_edge) has_to_be_changed(la)=.TRUE.
      END DO

      DO la=0, am_max
        IF(nf_per_l(la)==1) THEN
          ipos=ipos+1
          new_expo(ipos)=max_min_exp_per_l(1,la)
          IF(new_expo(ipos)>=exp_limits(1,ikind).AND.new_expo(ipos)<=exp_limits(2,ikind)) THEN
            max_rel_dev(ipos)=(new_expo(ipos)-old_lower_B(ipos))/new_expo(ipos)
            IF(max_rel_dev(ipos)<=0.1_dp) max_rel_dev(ipos)=0.8_dp
          ELSE
            new_expo(ipos)=(exp_limits(2,ikind)+exp_limits(1,ikind))/2.0_dp
            max_rel_dev(ipos)=(new_expo(ipos)-exp_limits(1,ikind))/new_expo(ipos)
          END IF
          IF(has_to_be_changed(la)) change_expo(ipos)=1
        ELSE
          IF(max_min_exp_per_l(2,la)/max_min_exp_per_l(1,la)<2.0_dp) THEN
            max_min_exp_per_l(1,la)=max_min_exp_per_l(1,la)*0.5
            max_min_exp_per_l(2,la)=max_min_exp_per_l(2,la)*1.5
          END IF
          geom_fact=(max_min_exp_per_l(2,la)/max_min_exp_per_l(1,la))**(1.0_dp/REAL(nf_per_l(la)-1,dp))
          DO iexpo=1, nf_per_l(la)
            ipos=ipos+1
            new_expo(ipos)=max_min_exp_per_l(1,la)*(geom_fact**(iexpo-1))
            max_rel_dev(ipos)=(geom_fact-1.0_dp)/(geom_fact+1.0_dp)*0.9_dp
            IF(has_to_be_changed(la)) change_expo(ipos)=1
          END DO
        END IF
      END DO

      DEALLOCATE(has_to_be_changed)

      DEALLOCATE(nf_per_l)
      DEALLOCATE(max_min_exp_per_l)
    END DO

    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        RI_basis_parameter(ikind)%zet(1,iset)=new_expo(ipos)
      END DO
    END DO

    CALL init_transf(nkind,RI_basis_parameter,lower_B,max_dev,max_rel_dev)

    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        IF(change_expo(ipos)==0) THEN
          ! restore original
          pnew(ipos) = old_pnew(ipos)
          lower_B(ipos) = old_lower_B(ipos)
          max_dev(ipos) = old_max_dev(ipos)
          max_rel_dev(ipos) = old_max_rel_dev(ipos)
          RI_basis_parameter(ikind)%zet(1,iset)=old_expo(ipos)
        END IF
      END DO
    END DO

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param DI ...
!> \param Emp2_RI ...
!> \param DRI ...
!> \param Emp2 ...
!> \param homo ...
!> \param homo_beta ...
!> \param virtual ...
!> \param virtual_beta ...
!> \param fact ...
!> \param fact2 ...
!> \param calc_ex ...
!> \param MOenerg ...
!> \param MOenerg_beta ...
!> \param abij ...
!> \param Lai ...
!> \param Lai_beta ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE contract_integrals(DI,Emp2_RI,DRI,Emp2,homo,homo_beta,virtual,virtual_beta,&
                                fact,fact2,calc_ex,&
                                MOenerg,MOenerg_beta,abij,&
                                Lai,Lai_beta,para_env,error)
    REAL(KIND=dp)                            :: DI, Emp2_RI, DRI, Emp2
    INTEGER                                  :: homo, homo_beta, virtual, &
                                                virtual_beta
    REAL(KIND=dp)                            :: fact, fact2
    LOGICAL                                  :: calc_ex
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MOenerg, MOenerg_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: abij
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Lai, Lai_beta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, b, i, ij_counter, j, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: t_iajb, t_iajb_RI
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mat_ab

    failure=.FALSE.

    ALLOCATE(mat_ab(virtual,virtual_beta),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
    DI=0.0_dp
    Emp2_RI=0.0_dp
    ij_counter=0
    DO j=1, homo_beta
      DO i=1, homo
        ij_counter=ij_counter+1
        IF (MOD(ij_counter,para_env%num_pe)/=para_env%mepos) CYCLE
        mat_ab=0.0_dp
        mat_ab(:,:)=MATMUL(TRANSPOSE(Lai(:,:,i)),Lai_beta(:,:,j))
        DO b=1, virtual_beta
          DO a=1, virtual
            IF(calc_ex) THEN
              t_iajb=fact*abij(a,b,i,j)-abij(b,a,i,j)
              t_iajb_RI=fact*mat_ab(a,b)-mat_ab(b,a)
            ELSE
              t_iajb=fact*abij(a,b,i,j)
              t_iajb_RI=fact*mat_ab(a,b)
            END IF
            t_iajb=t_iajb/(MOenerg(a+homo)+MOenerg_beta(b+homo_beta)-MOenerg(i)-MOenerg_beta(j))
            t_iajb_RI=t_iajb_RI/(MOenerg(a+homo)+MOenerg_beta(b+homo_beta)-MOenerg(i)-MOenerg_beta(j))

            Emp2_RI=Emp2_RI-t_iajb_RI*mat_ab(a,b)*fact2

            DI=DI-t_iajb*mat_ab(a,b)*fact2

          END DO
        END DO
      END DO
    END DO
    CALL mp_sum(DI,para_env%group)
    CALL mp_sum(Emp2_RI,para_env%group)

    DRI=Emp2-Emp2_RI
    DI=2.0D+00*DI-Emp2-Emp2_RI

    DEALLOCATE(mat_ab)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param homo ...
!> \param homo_beta ...
!> \param para_env ...
!> \param elements_ij_proc ...
!> \param ij_list_proc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calc_elem_ij_proc(homo,homo_beta,para_env,elements_ij_proc,ij_list_proc,error)
    INTEGER                                  :: homo, homo_beta
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: elements_ij_proc
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_list_proc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_elem_ij_proc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij_counter, j, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    elements_ij_proc=0
    ij_counter=-1
    DO i=1, homo
      DO j=1, homo_beta
        ij_counter=ij_counter+1
        IF (MOD(ij_counter,para_env%num_pe)==para_env%mepos) elements_ij_proc=elements_ij_proc+1
      END DO
    END DO

    ALLOCATE(ij_list_proc(elements_ij_proc,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ij_list_proc=0
    ij_counter=-1
    elements_ij_proc=0
    DO i=1, homo
      DO j=1, homo_beta
        ij_counter=ij_counter+1
        IF (MOD(ij_counter,para_env%num_pe)==para_env%mepos) THEN
           elements_ij_proc=elements_ij_proc+1
           ij_list_proc(elements_ij_proc,1)=i
           ij_list_proc(elements_ij_proc,2)=j
        END IF
      END DO
    END DO

  END SUBROUTINE calc_elem_ij_proc

! *****************************************************************************
!> \brief ...
!> \param lower_B ...
!> \param max_dev ...
!> \param valin ...
!> \param valout ...
! *****************************************************************************
  SUBROUTINE transf_val(lower_B, max_dev, valin, valout)
    REAL(KIND=dp)                            :: lower_B, max_dev, valin, &
                                                valout

    REAL(KIND=dp), PARAMETER                 :: alpha = 2.633915794_dp

    valout=0.0_dp
    valout=lower_B+max_dev/(1.0_dp+EXP(-alpha*valin))

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param nkind ...
!> \param max_rel_dev_output ...
!> \param basis_was_assoc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE generate_RI_init_basis(qs_env,mp2_env,nkind,max_rel_dev_output,basis_was_assoc,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    INTEGER                                  :: nkind
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: max_rel_dev_output
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: basis_was_assoc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_RI_init_basis', &
      routineP = moduleN//':'//routineN

    INTEGER :: basis_quality, handle, iexpo, iii, ikind, ipgf, iset, ishell, &
      jexpo, jjj, la, max_am, nexpo_shell, nseta, nsgfa, nsgfa_RI, prog_func, &
      prog_l, RI_max_am, RI_nset, RI_prev_size, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: l_expo, num_sgf_per_l, &
                                                ordered_pos, RI_l_expo, &
                                                RI_num_sgf_per_l, &
                                                tot_num_exp_per_l
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: l_tab
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, nl
    LOGICAL                                  :: external_num_of_func, failure
    REAL(dp)                                 :: geom_fact
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: exponents, RI_exponents
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: exp_tab, max_min_exp_l
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_a, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: gcc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: max_rel_dev, max_rel_dev_prev
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, tmp_basis
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: atom_kind

    CALL timeset(routineN,handle)
    failure = .FALSE.

    basis_quality=mp2_env%ri_opt_param%basis_quality
    external_num_of_func=.FALSE.
    IF(ALLOCATED(mp2_env%ri_opt_param%RI_nset_per_l)) external_num_of_func=.TRUE.

    NULLIFY(qs_kind_set)
    CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, error=error)
    nkind = SIZE(qs_kind_set,1) 

    ALLOCATE(basis_was_assoc(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    basis_was_assoc=.FALSE.

    IF(external_num_of_func.AND.nkind>1) THEN
      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                     "There are more than one kind of atom. The same pattern of functions, "//&
                     "as specified by NUM_FUNC, will be used for all kinds.", &
                     error,only_ionode=.TRUE.)
    END IF

    DO ikind = 1, nkind
      NULLIFY(atom_kind)
      atom_kind => qs_kind_set(ikind)

      CALL get_qs_kind(qs_kind=atom_kind,basis_set=orb_basis_a,basis_type="RI_AUX",error=error)
      ! save info if the basis was or not associated
      basis_was_assoc(ikind)=ASSOCIATED(orb_basis_a)

      ! skip the generation of the initial guess if the RI basis is 
      ! provided in input
      IF(basis_was_assoc(ikind)) THEN
        nseta = orb_basis_a%nset
        la_max => orb_basis_a%lmax
        la_min => orb_basis_a%lmin
        npgfa => orb_basis_a%npgf
        nshell => orb_basis_a%nshell
        zet => orb_basis_a%zet
        nl => orb_basis_a%l

        max_am=MAXVAL(la_max)

        RI_max_am=max_am
        ALLOCATE(RI_num_sgf_per_l(0:RI_max_am),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        RI_num_sgf_per_l=0
        RI_nset=0
        DO iset=1, nseta
          DO la=la_min(iset), la_max(iset)
            RI_nset=RI_nset+1
            RI_num_sgf_per_l(la)=RI_num_sgf_per_l(la)+1
            IF(npgfa(iset)>1) THEN
              CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                   "The RI basis set optimizer can not handle contracted Gaussian. "//&
                   "Calculation continue with only uncontracted functions.", &
                   error,only_ionode=.TRUE.)
            END IF
          END DO
        END DO

        ALLOCATE(exp_tab(MAXVAL(RI_num_sgf_per_l),0:RI_max_am),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        exp_tab=0.0_dp
        DO iii=0, RI_max_am
          iexpo=0
          DO iset=1, nseta
            DO la=la_min(iset), la_max(iset)
              IF(la/=iii) CYCLE
              iexpo=iexpo+1
              exp_tab(iexpo,iii)=zet(1,iset) 
            END DO
          END DO
        END DO

        ! sort exponents
        DO iii=0, RI_max_am
          ALLOCATE(ordered_pos(RI_num_sgf_per_l(iii)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ordered_pos=0
          CALL sort(exp_tab(1:RI_num_sgf_per_l(iii),iii), RI_num_sgf_per_l(iii), ordered_pos)
          DEALLOCATE(ordered_pos)
        END DO

        ALLOCATE(RI_l_expo(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        RI_l_expo=0
        ALLOCATE(RI_exponents(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        RI_exponents=0.0_dp

        iset=0
        DO iii=0, RI_max_am
          DO iexpo=1, RI_num_sgf_per_l(iii)
            iset=iset+1
            RI_l_expo(iset)=iii
            RI_exponents(iset)=exp_tab(iexpo,iii)
          END DO
        END DO
        DEALLOCATE(exp_tab)


        ALLOCATE(max_rel_dev(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        max_rel_dev=1.0_dp
        iset=0
        DO iii=0, RI_max_am
          IF(RI_num_sgf_per_l(iii)==1) THEN
            iset=iset+1
            max_rel_dev(iset)=0.35_dp
          ELSE
            iset=iset+1
            max_rel_dev(iset)=(RI_exponents(iset+1)+RI_exponents(iset))/2.0_dp
            max_rel_dev(iset)=max_rel_dev(iset)/RI_exponents(iset)-1.0_dp
            DO iexpo=2, RI_num_sgf_per_l(iii)
              iset=iset+1
              max_rel_dev(iset)=(RI_exponents(iset)+RI_exponents(iset-1))/2.0_dp
              max_rel_dev(iset)=1.0_dp-max_rel_dev(iset)/RI_exponents(iset)
            END DO
          END IF
        END DO
        max_rel_dev(:) = max_rel_dev(:)*0.9_dp
        DEALLOCATE(RI_num_sgf_per_l)

        ! deallocate the old basis before moving on
        CALL remove_basis_from_container(qs_kind_set(ikind)%basis_sets,basis_type="RI_AUX",error=error)

      ELSE
 
        CALL get_qs_kind(qs_kind=atom_kind,basis_set=orb_basis_a,error=error)

        sphi_a => orb_basis_a%sphi
        nseta = orb_basis_a%nset
        la_max => orb_basis_a%lmax
        la_min => orb_basis_a%lmin
        npgfa => orb_basis_a%npgf
        first_sgfa => orb_basis_a%first_sgf
        nshell => orb_basis_a%nshell
        zet => orb_basis_a%zet
        gcc => orb_basis_a%gcc
        nl => orb_basis_a%l
        nsgfa = orb_basis_a%nsgf
        
        max_am=MAXVAL(la_max)

        nexpo_shell=0
        DO iset = 1, nseta
          DO ishell=1,nshell(iset)
            DO la = la_min(iset), la_max(iset)
              IF(ishell>1) THEN
                IF(nl(ishell,iset)==nl(ishell-1,iset)) CYCLE
              END IF
              IF(la/=nl(ishell,iset)) CYCLE
              nexpo_shell=nexpo_shell+npgfa(iset)
            END DO
          END DO
        END DO  

        ALLOCATE(exponents(nexpo_shell),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        exponents=0.0_dp
        ALLOCATE(l_expo(nexpo_shell),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        l_expo=0
        ALLOCATE(num_sgf_per_l(0:max_am),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        num_sgf_per_l=0
        iexpo=0
        DO iset = 1, nseta
          DO ishell=1, nshell(iset)
            DO la = la_min(iset), la_max(iset)
              IF(ishell>1) THEN
                IF(nl(ishell,iset)==nl(ishell-1,iset)) CYCLE
              END IF
              IF(la/=nl(ishell,iset)) CYCLE
              DO ipgf=1, npgfa(iset)
                iexpo=iexpo+1
                exponents(iexpo)=zet(ipgf,iset)
                l_expo(iexpo)=la
              END DO
              num_sgf_per_l(la)=num_sgf_per_l(la)+1
            END DO
          END DO
        END DO

        ALLOCATE(exp_tab(nexpo_shell,nexpo_shell),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        exp_tab=0.0_dp
        ALLOCATE(l_tab(nexpo_shell,nexpo_shell),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        l_tab=0
        ALLOCATE(tot_num_exp_per_l(0:max_am*2),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        tot_num_exp_per_l=0
        DO iexpo=1, nexpo_shell
          DO jexpo=iexpo, nexpo_shell
            exp_tab(jexpo,iexpo)=exponents(jexpo)+exponents(iexpo)
            exp_tab(iexpo,jexpo)=exp_tab(jexpo,iexpo)
            l_tab(jexpo,iexpo)=l_expo(jexpo)+l_expo(iexpo)
            l_tab(iexpo,jexpo)=l_tab(jexpo,iexpo)
            tot_num_exp_per_l(l_tab(jexpo,iexpo))=tot_num_exp_per_l(l_tab(jexpo,iexpo))+1
          END DO
        END DO
        DEALLOCATE(l_expo)
        DEALLOCATE(exponents)

        ALLOCATE(max_min_exp_l(2,0:max_am*2),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        max_min_exp_l(1,:)=HUGE(0)
        max_min_exp_l(2,:)=-HUGE(0)

        DO la=0, max_am*2
          DO iexpo=1, nexpo_shell
            DO jexpo=iexpo, nexpo_shell
              IF(la/=l_tab(jexpo,iexpo)) CYCLE
              IF(exp_tab(jexpo,iexpo)<=max_min_exp_l(1,la)) max_min_exp_l(1,la)=exp_tab(jexpo,iexpo)
              IF(exp_tab(jexpo,iexpo)>=max_min_exp_l(2,la)) max_min_exp_l(2,la)=exp_tab(jexpo,iexpo)
            END DO
          END DO
        END DO
        DEALLOCATE(exp_tab)
        DEALLOCATE(l_tab)
        
        ! ! scale the limits of the exponents to avoid the optimizer to go out-of-range 
        ! ! (0.2 just empirical parameter)
        max_min_exp_l(1,:)=max_min_exp_l(1,:)/0.80_dp
        max_min_exp_l(2,:)=max_min_exp_l(2,:)/1.20_dp

        ALLOCATE(RI_num_sgf_per_l(0:max_am*2),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        RI_num_sgf_per_l=0
        
        SELECT CASE(basis_quality)
          CASE(0) 
            ! normal
            prog_func=0
            prog_l=2
          CASE(1)
            ! large
            prog_func=1
            prog_l=2
          CASE(2)
            prog_func=2
            prog_l=2
          CASE DEFAULT
            prog_func=0
            prog_l=2
        END SELECT

        IF(external_num_of_func) THEN
          ! cp2k can not exceed angular momentum 7
          RI_max_am=MIN(SIZE(mp2_env%ri_opt_param%RI_nset_per_l)-1,7)
          IF(RI_max_am>max_am*2) THEN
            DEALLOCATE(RI_num_sgf_per_l)
            ALLOCATE(RI_num_sgf_per_l(0:RI_max_am),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            RI_num_sgf_per_l=0
          END IF
          DO la=0, RI_max_am
            RI_num_sgf_per_l(la)=mp2_env%ri_opt_param%RI_nset_per_l(la)
          END DO
        ELSE
          RI_num_sgf_per_l(0)=num_sgf_per_l(0)*2+prog_func
          DO la=1, max_am
            RI_num_sgf_per_l(la)=RI_num_sgf_per_l(la-1)-1
            IF(RI_num_sgf_per_l(la)==0) THEN
              RI_num_sgf_per_l(la)=1
              EXIT
            END IF
          END DO
          DO la=max_am+1, max_am*2
            RI_num_sgf_per_l(la)=RI_num_sgf_per_l(la-1)-prog_l
            IF(RI_num_sgf_per_l(la)==0) THEN
              RI_num_sgf_per_l(la)=1 
            END IF
            IF(RI_num_sgf_per_l(la)==1) EXIT
          END DO
          RI_max_am=MIN(max_am*2,7)
          DO la=0, MIN(max_am*2,7)
            IF(RI_num_sgf_per_l(la)==0) THEN
              RI_max_am=la-1 
              EXIT
            END IF
          END DO

          iii=0
          jjj=0
          nsgfa_RI=0
          DO la=1, max_am*2
            IF(tot_num_exp_per_l(la)>=iii) THEN
              iii=tot_num_exp_per_l(la)
              jjj=la
            END IF
            nsgfa_RI=nsgfa_RI+RI_num_sgf_per_l(la)*(la*2+1)
          END DO
          DEALLOCATE(tot_num_exp_per_l)
          IF(REAL(nsgfa_RI,KIND=dp)/REAL(nsgfa,KIND=dp)<=2.5_dp) THEN
            RI_num_sgf_per_l(jjj)=RI_num_sgf_per_l(jjj)+1
          END IF
        END IF

        RI_nset=SUM(RI_num_sgf_per_l)

        ALLOCATE(RI_exponents(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        RI_exponents=0.0_dp

        ALLOCATE(RI_l_expo(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        RI_l_expo=0

        ALLOCATE(max_rel_dev(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        max_rel_dev=1.0_dp

        iset=0
        DO la=0, RI_max_am
          IF(RI_num_sgf_per_l(la)==1) THEN
            iset=iset+1
            RI_exponents(iset)=(max_min_exp_l(2,la)+max_min_exp_l(1,la))/2.0_dp
            RI_l_expo(iset)=la
            geom_fact=max_min_exp_l(2,la)/max_min_exp_l(1,la)

            max_rel_dev(iset)=(geom_fact-1.0_dp)/(geom_fact+1.0_dp)*0.9_dp
          ELSE
            geom_fact=(max_min_exp_l(2,la)/max_min_exp_l(1,la))**(1.0_dp/REAL(RI_num_sgf_per_l(la)-1,dp))
            DO iexpo=1, RI_num_sgf_per_l(la)
              iset=iset+1
              RI_exponents(iset)=max_min_exp_l(1,la)*(geom_fact**(iexpo-1))
              RI_l_expo(iset)=la
              max_rel_dev(iset)=(geom_fact-1.0_dp)/(geom_fact+1.0_dp)*0.9_dp
            END DO
          END IF
        END DO
        DEALLOCATE(num_sgf_per_l)
        DEALLOCATE(max_min_exp_l)
        DEALLOCATE(RI_num_sgf_per_l)

        ! IF(qs_env%para_env%mepos==0) THEN
        !   WRITE(*,*) RI_nset, RI_max_am, geom_fact
        !   iset=0
        !   DO la=0, RI_max_am
        !     DO iexpo=1, RI_num_sgf_per_l(la)
        !       iset=iset+1
        !       WRITE(*,*) RI_l_expo(iset), RI_exponents(iset)
        !     END DO
        !   END DO
        !   WRITE(*,*) nsgfa, nsgfa_RI
        ! END IF
 
      END IF ! RI basis not associated

      ! create the new basis
      NULLIFY(tmp_basis)
      CALL create_ri_basis(tmp_basis,RI_nset,RI_l_expo,RI_exponents,error)
      CALL add_basis_set_to_container(qs_kind_set(ikind)%basis_sets,tmp_basis,"RI_AUX",error=error)
!d    CALL copy_gto_basis_set(tmp_basis,qs_kind_set(ikind)%ri_aux_basis_set,error)

      DEALLOCATE(RI_exponents)
      DEALLOCATE(RI_l_expo)

      IF(.NOT.ALLOCATED(max_rel_dev_output)) THEN
        ALLOCATE(max_rel_dev_output(RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        max_rel_dev_output(:)=max_rel_dev
      ELSE
        ! make a copy
        RI_prev_size=SIZE(max_rel_dev_output)
        ALLOCATE(max_rel_dev_prev(RI_prev_size),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        max_rel_dev_prev(:)=max_rel_dev_output
        DEALLOCATE(max_rel_dev_output)
        ! reallocate and copy
        ALLOCATE(max_rel_dev_output(RI_prev_size+RI_nset),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        max_rel_dev_output(1:RI_prev_size)=max_rel_dev_prev
        max_rel_dev_output(RI_prev_size+1:RI_prev_size+RI_nset)=max_rel_dev
        DEALLOCATE(max_rel_dev_prev)
      END IF
      DEALLOCATE(max_rel_dev)

    END DO ! ikind

    IF(external_num_of_func) DEALLOCATE(mp2_env%ri_opt_param%RI_nset_per_l)

    CALL timestop(handle)   

  END SUBROUTINE generate_RI_init_basis

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
!> \param RI_nset ...
!> \param RI_l_expo ...
!> \param RI_exponents ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_ri_basis(gto_basis_set,RI_nset,RI_l_expo,RI_exponents,error)
    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    INTEGER                                  :: RI_nset
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: RI_l_expo
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: RI_exponents
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_ri_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ico, ipgf, iset, &
                                                ishell, lshell, m, maxco, &
                                                maxl, maxpgf, maxshell, ncgf, &
                                                nmin, nset, nsgf, stat
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, n
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(lmax, lmin, npgf, nshell, l, n, zet, gcc)

    ! allocate the basis
    CALL allocate_gto_basis_set(gto_basis_set,error)

    ! brute force
    nset = 0
    maxshell=0
    maxpgf=0
    maxco=0
    ncgf=0
    nsgf=0
    gto_basis_set%nset = nset
    gto_basis_set%ncgf = ncgf
    gto_basis_set%nsgf = nsgf
    CALL reallocate(gto_basis_set%lmax,1,nset)
    CALL reallocate(gto_basis_set%lmin,1,nset)
    CALL reallocate(gto_basis_set%npgf,1,nset)
    CALL reallocate(gto_basis_set%nshell,1,nset)
    CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%set_radius,1,nset)
    CALL reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%ncgf_set,1,nset)
    CALL reallocate(gto_basis_set%nsgf_set,1,nset)
    CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
    CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%lx,1,ncgf)
    CALL reallocate(gto_basis_set%ly,1,ncgf)
    CALL reallocate(gto_basis_set%lz,1,ncgf)
    CALL reallocate(gto_basis_set%m,1,nsgf)
    CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)

    nset=RI_nset

    ! locals
    CALL reallocate(npgf,1,nset)
    CALL reallocate(nshell,1,nset)
    CALL reallocate(lmax,1,nset)
    CALL reallocate(lmin,1,nset)
    CALL reallocate(n,1,1,1,nset)

    maxl = 0
    maxpgf = 0
    maxshell = 0
    DO iset=1,nset
      n(1,iset)=iset
      lmin(iset)=RI_l_expo(iset)
      lmax(iset)=RI_l_expo(iset)
      npgf(iset)=1

      maxl = MAX(maxl,lmax(iset))

      IF (npgf(iset) > maxpgf) THEN
         maxpgf = npgf(iset)
         CALL reallocate(zet,1,maxpgf,1,nset)
         CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
      END IF
      nshell(iset) = 0
      DO lshell=lmin(iset),lmax(iset)
         nmin = n(1,iset) + lshell - lmin(iset)
         ishell=1
         nshell(iset) = nshell(iset) + ishell
         IF (nshell(iset) > maxshell) THEN
            maxshell = nshell(iset)
            CALL reallocate(n,1,maxshell,1,nset)
            CALL reallocate(l,1,maxshell,1,nset)
            CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
         END IF
         DO i=1,ishell
            n(nshell(iset)-ishell+i,iset) = nmin + i - 1
            l(nshell(iset)-ishell+i,iset) = lshell
         END DO
      END DO

      DO ipgf=1,npgf(iset)
        zet(ipgf,iset)=RI_exponents(iset)
        DO ishell=1, nshell(iset)
          gcc(ipgf,ishell,iset)=1.0_dp
        END DO
      END DO
    END DO

    CALL init_orbital_pointers(maxl)

    gto_basis_set%nset = nset
    CALL reallocate(gto_basis_set%lmax,1,nset)
    CALL reallocate(gto_basis_set%lmin,1,nset)
    CALL reallocate(gto_basis_set%npgf,1,nset)
    CALL reallocate(gto_basis_set%nshell,1,nset)
    CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)

    ! Copy the basis set information into the data structure

    DO iset=1,nset
       gto_basis_set%lmax(iset) = lmax(iset)
       gto_basis_set%lmin(iset) = lmin(iset)
       gto_basis_set%npgf(iset) = npgf(iset)
       gto_basis_set%nshell(iset) = nshell(iset)
       DO ishell=1,nshell(iset)
          gto_basis_set%n(ishell,iset) = n(ishell,iset)
          gto_basis_set%l(ishell,iset) = l(ishell,iset)
          DO ipgf=1,npgf(iset)
             gto_basis_set%gcc(ipgf,ishell,iset) = gcc(ipgf,ishell,iset)
          END DO
       END DO
       DO ipgf=1,npgf(iset)
          gto_basis_set%zet(ipgf,iset) = zet(ipgf,iset)
       END DO
    END DO

    ! Initialise the depending atomic kind information

    CALL reallocate(gto_basis_set%set_radius,1,nset)
    CALL reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%ncgf_set,1,nset)
    CALL reallocate(gto_basis_set%nsgf_set,1,nset)

    maxco = 0
    ncgf = 0
    nsgf = 0

    DO iset=1,nset
       gto_basis_set%ncgf_set(iset) = 0
       gto_basis_set%nsgf_set(iset) = 0
       DO ishell=1,nshell(iset)
          lshell = gto_basis_set%l(ishell,iset)
          gto_basis_set%first_cgf(ishell,iset) = ncgf + 1
          ncgf = ncgf + nco(lshell)
          gto_basis_set%last_cgf(ishell,iset) = ncgf
          gto_basis_set%ncgf_set(iset) =&
               gto_basis_set%ncgf_set(iset) + nco(lshell)
          gto_basis_set%first_sgf(ishell,iset) = nsgf + 1
          nsgf = nsgf + nso(lshell)
          gto_basis_set%last_sgf(ishell,iset) = nsgf
          gto_basis_set%nsgf_set(iset) =&
               gto_basis_set%nsgf_set(iset) + nso(lshell)
       END DO
       maxco = MAX(maxco,npgf(iset)*ncoset(lmax(iset)))
    END DO

    gto_basis_set%ncgf = ncgf
    gto_basis_set%nsgf = nsgf

    CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
    CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%lx,1,ncgf)
    CALL reallocate(gto_basis_set%ly,1,ncgf)
    CALL reallocate(gto_basis_set%lz,1,ncgf)
    CALL reallocate(gto_basis_set%m,1,nsgf)
    CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)
    ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ncgf = 0
    nsgf = 0

    DO iset=1,nset
       DO ishell=1,nshell(iset)
          lshell = gto_basis_set%l(ishell,iset)
          DO ico=ncoset(lshell-1)+1,ncoset(lshell)
             ncgf = ncgf + 1
             gto_basis_set%lx(ncgf) = indco(1,ico)
             gto_basis_set%ly(ncgf) = indco(2,ico)
             gto_basis_set%lz(ncgf) = indco(3,ico)
             gto_basis_set%cgf_symbol(ncgf) =&
                  cgf_symbol(n(ishell,iset),(/gto_basis_set%lx(ncgf),&
                  gto_basis_set%ly(ncgf),&
                  gto_basis_set%lz(ncgf)/))
          END DO
          DO m=-lshell,lshell
             nsgf = nsgf + 1
             gto_basis_set%m(nsgf) = m
             gto_basis_set%sgf_symbol(nsgf) =&
                  sgf_symbol(n(ishell,iset),lshell,m)
          END DO
       END DO
    END DO

    DEALLOCATE (gcc,l,lmax,lmin,n,npgf,nshell,zet,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE

END MODULE mp2_optimize_ri_basis

