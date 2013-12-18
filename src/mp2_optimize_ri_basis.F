!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
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
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
! USE f77_blas
  USE hfx_types,                       ONLY: hfx_basis_info_type,&
                                             hfx_basis_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_comm_split_direct,&
                                             mp_sum
  USE mp2_direct_method,               ONLY: mp2_canonical_direct_single_batch
  USE mp2_ri_libint,                   ONLY: libint_ri_mp2,&
                                             read_RI_basis_set,&
                                             release_RI_basis_set
  USE mp2_types,                       ONLY: mp2_biel_type,&
                                             mp2_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_optimize_ri_basis'

  PUBLIC :: optimize_ri_basis_main

  CONTAINS

! *****************************************************************************
!> \brief optimize RI-MP2 basis set
!> \author Mauro Del Ben
! *****************************************************************************
  SUBROUTINE optimize_ri_basis_main(Emp2,Emp2_Cou,Emp2_ex,Emp2_S,Emp2_T,dimen,natom,homo, &
                                     mp2_biel,mp2_env,C,Auto, &
                                     kind_of,basis_parameter, &
                                     qs_env,particle_set,matrix_ks,rho,hfx_sections,para_env, &
                                     unit_nr,error,&
                                     homo_beta,C_beta,Auto_beta)

    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_ex, &
                                                Emp2_S, Emp2_T
    INTEGER                                  :: dimen, natom, homo
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
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
      iiter, ikind, ipgf, iset, ishell, j, local_unit_nr, max_num_iter, ndof, &
      nkind, number_groups, stat, virtual, virtual_beta
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_list_proc, index_table_RI
    LOGICAL                                  :: failure, open_shell_case
    REAL(KIND=dp) :: DI, DI_new, DRI, DRI_new, Emp2_AA, Emp2_AA_Cou, &
      Emp2_AA_ex, Emp2_AB, Emp2_AB_Cou, Emp2_AB_ex, Emp2_BB, Emp2_BB_Cou, &
      Emp2_BB_ex, Emp2_RI, Emp2_RI_new, eps_DI_rel, eps_DRI, eps_step, fac, &
      fad, fae, sumdg, sumxi
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: deriv, dg, g, hdg, p, pnew, xi
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: hessin
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: Integ_MP2, Integ_MP2_AA, &
                                                Integ_MP2_AB, Integ_MP2_BB
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type)                      :: error_sub
    TYPE(cp_logger_type), POINTER            :: logger, logger_sub
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    TYPE(gto_basis_set_type), DIMENSION(:), &
      POINTER                                :: basis_set_RI
    TYPE(hfx_basis_info_type)                :: RI_basis_info
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_S0, RI_basis_parameter

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
      CALL calc_elem_ij_proc(homo,homo,para_env,elements_ij_proc,ij_list_proc)
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
      CALL calc_elem_ij_proc(homo_beta,homo_beta,para_env,elements_ij_proc,ij_list_proc)
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
      CALL calc_elem_ij_proc(homo,homo_beta,para_env,elements_ij_proc,ij_list_proc)
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
      CALL calc_elem_ij_proc(homo,homo,para_env,elements_ij_proc,ij_list_proc)
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
    CALL cp_error_init(error_sub, stop_level=cp_failure_level, logger=logger_sub)

    CALL read_RI_basis_set(qs_env,RI_basis_parameter,RI_basis_info,&
                           natom,nkind,kind_of,index_table_RI,dimen_RI,&
                           basis_S0,error) 
    ndof=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ndof=ndof+1
      END DO
    END DO

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
    ALLOCATE(hessin(ndof,ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    hessin=0.0_dp
    DO i=1, ndof 
      hessin(i,i)=1.0_dp
    END DO

    ALLOCATE(deriv(ndof),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    deriv=0.0_dp

    ! Get the position vector
    CALL basis2p(nkind,ndof,basis_set_RI,RI_basis_parameter,p)

    ! Calculate RI-MO-ERI's
    CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                          Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                          qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                          kind_of,index_table_RI,basis_set_RI,mp2_biel,mp2_env,Auto,C,&
                          hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                          open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                          .TRUE.,error)

    ! ! Calculate function (DI) derivatives with respect to the RI basis exponent
    CALL calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps_step,&
                              qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,basis_set_RI,mp2_biel,mp2_env,Auto,C,&
                              hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                              para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                              deriv,error,error_sub)
    g=deriv
    xi=-g

    ! get the atomic kind set for writing the basis
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)

    DO iiter=1, max_num_iter
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,I5)') 'OPTIMIZATION STEP NUMBER',iiter

      ! perform step
      pnew=p+xi
      CALL p2basis(nkind,ndof,basis_set_RI,RI_basis_parameter,pnew)

      ! calculate energy at the new point
      CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI_new,DRI_new,DI_new,&
                            Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                            qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                            kind_of,index_table_RI,basis_set_RI,mp2_biel,mp2_env,Auto,C,&
                            hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                            open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                            .FALSE.,error)
 
      ! update energy and direction
      DI=DI_new
      xi=pnew-p
      p=pnew

      ! check for convergence
      IF (unit_nr>0) THEN
        WRITE(unit_nr,*)
        DO ikind=1, nkind
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
                                             (atomic_kind_set(ikind)%ri_aux_basis_set%&
                                             gcc(ipgf,ishell,iset),&
                                             ishell=1,atomic_kind_set(ikind)%ri_aux_basis_set%nshell(iset))
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
      CALL calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                                Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps_step,&
                                qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                                kind_of,index_table_RI,basis_set_RI,mp2_biel,mp2_env,Auto,C,&
                                hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                                open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                                para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                                deriv,error,error_sub)

      ! g is the vector containing the old gradient
      dg=deriv-g
      g=deriv
      hdg=MATMUL(hessin,dg)

      fac=SUM(dg*xi)
      fae=SUM(dg*hdg)
      sumdg=SUM(dg*dg)
      sumxi=SUM(xi*xi)
      
      IF(fac**2>sumdg*sumxi*3.0E-8_dp) THEN
        fac=1.0_dp/fac
        fad=1.0_dp/fae
        dg=fac*xi-fad*hdg
        DO i=1, ndof
          DO j=1, ndof
            hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)&
                                   -fad*hdg(i)*hdg(j)&
                                   +fae*dg(i)*dg(j)
          END DO
        END DO
      ELSE
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Skip Hessian Update'   
      END IF

      !XXXXXXXXXXXXX
      ! DO i=1, ndof
      !   IF(g(i)<0.0D+00) THEN
      !     xi(i)=pnew(i)*0.005
      !   ELSE
      !     xi(i)=-pnew(i)*0.005
      !   END IF
      ! END DO
      !XXXXXXXXXXXXX

      ! new direction
      xi=-MATMUL(hessin,g)

    END DO

    IF(.NOT.(DI/ABS(Emp2)<=eps_DI_rel.AND.ABS(DRI_new)<=eps_DRI)) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,I5,A)') 'OPTIMIZATION NOT CONVERGED IN',max_num_iter,' STEPS.'
      IF (unit_nr>0) WRITE(unit_nr,*)
    END IF

    DEALLOCATE(p)
    DEALLOCATE(xi)
    DEALLOCATE(g)
    DEALLOCATE(pnew)
    DEALLOCATE(dg)  
    DEALLOCATE(hdg)  
    DEALLOCATE(Hessin)

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

  SUBROUTINE calc_energy_func_der(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI_ref,DRI_ref,DI_ref,&
                                  Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,eps,&
                                  qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                                  kind_of,index_table_RI,basis_set_RI_original,mp2_biel,mp2_env,Auto,C,&
                                  hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                                  open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                                  para_env,para_env_sub,number_groups,color_sub,unit_nr,&
                                  deriv,error,error_sub)
    REAL(KIND=dp)                            :: Emp2, Emp2_AA, Emp2_BB, &
                                                Emp2_AB, Emp2_RI_ref, &
                                                DRI_ref, DI_ref
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: Integ_MP2, Integ_MP2_AA, &
                                                Integ_MP2_BB, Integ_MP2_AB
    REAL(KIND=dp)                            :: eps
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: nkind, natom, dimen, &
                                                dimen_RI, homo, virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: index_table_RI
    TYPE(gto_basis_set_type), DIMENSION(:), &
      POINTER                                :: basis_set_RI_original
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter, &
                                                RI_basis_parameter
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
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: deriv
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_energy_func_der', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ideriv, ikind, iset, &
                                                nseta
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: DI, DRI, Emp2_RI
    REAL(KIND=dp), VOLATILE                  :: step, temp
    TYPE(gto_basis_set_type), DIMENSION(:), &
      POINTER                                :: basis_set_RI

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! cycle over the RI basis set exponent
    deriv=0.0_dp
    ideriv=0
    DO ikind=1, nkind
      nseta=RI_basis_parameter(ikind)%nset
      DO iset=1, nseta
        ! for now only uncontracted aux basis set
        ideriv=ideriv+1
        IF (MOD(ideriv,number_groups)/=color_sub) CYCLE

        ! calculate the numerical derivative
        ! The eps is the relative change of the exponent for the
        ! calculation of the numerical derivative
        CPPostcondition(RI_basis_parameter(ikind)%npgf(iset)==1,cp_failure_level,routineP,error,failure)
        step=eps*RI_basis_parameter(ikind)%zet(1,iset)
        temp=RI_basis_parameter(ikind)%zet(1,iset)+step
        step=temp-RI_basis_parameter(ikind)%zet(1,iset)
        RI_basis_parameter(ikind)%zet(1,iset)=RI_basis_parameter(ikind)%zet(1,iset)+step
     
        CALL calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                              qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,basis_set_RI,mp2_biel,mp2_env,Auto,C,&
                              hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,&
                              para_env_sub,unit_nr,.TRUE.,error_sub)

        RI_basis_parameter(ikind)%zet(1,iset)=RI_basis_parameter(ikind)%zet(1,iset)-step

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

  SUBROUTINE calc_energy_func(Emp2,Emp2_AA,Emp2_BB,Emp2_AB,Emp2_RI,DRI,DI,&
                              Integ_MP2,Integ_MP2_AA,Integ_MP2_BB,Integ_MP2_AB,&
                              qs_env,particle_set,nkind,natom,dimen,dimen_RI,homo,virtual,&
                              kind_of,index_table_RI,basis_set_RI_original,mp2_biel,mp2_env,Auto,C,&
                              hfx_sections,basis_parameter,RI_basis_parameter,RI_basis_info,basis_S0,&
                              open_shell_case,homo_beta,virtual_beta,Auto_beta,C_beta,para_env,unit_nr,&
                              no_write,error)
    REAL(KIND=dp)                            :: Emp2, Emp2_AA, Emp2_BB, &
                                                Emp2_AB, Emp2_RI, DRI, DI
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: Integ_MP2, Integ_MP2_AA, &
                                                Integ_MP2_BB, Integ_MP2_AB
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: nkind, natom, dimen, &
                                                dimen_RI, homo, virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: index_table_RI
    TYPE(gto_basis_set_type), DIMENSION(:), &
      POINTER                                :: basis_set_RI_original
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter, &
                                                RI_basis_parameter
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
                       kind_of,basis_parameter,particle_set,&
                       RI_basis_parameter,RI_basis_info,basis_S0,index_table_RI,&
                       qs_env,hfx_sections,para_env,&
                       unit_nr,Lai,error)
    IF(open_shell_case) THEN
      CALL libint_ri_mp2(dimen,dimen_RI,homo_beta,natom,mp2_biel,mp2_env,C_beta,&
                         kind_of,basis_parameter,particle_set,&
                         RI_basis_parameter,RI_basis_info,basis_S0,index_table_RI,&
                         qs_env,hfx_sections,para_env,&
                         unit_nr,Lai_beta,error)
    END IF

    ! Contract integrals into energy
    IF(open_shell_case) THEN
      ! alpha-alpha
      CALL contract_integrals(DI_AA,Emp2_RI_AA,DRI_AA,Emp2_AA,dimen_RI,homo,homo,virtual,virtual,&
                              1.0_dp,0.5_dp,.TRUE.,&
                              Auto,Auto,Integ_MP2_AA,&
                              Lai,Lai,para_env)

      ! beta-beta
      CALL contract_integrals(DI_BB,Emp2_RI_BB,DRI_BB,Emp2_BB,dimen_RI,homo_beta,homo_beta,virtual_beta,virtual_beta,&
                              1.0_dp,0.5_dp,.TRUE.,&
                              Auto_beta,Auto_beta,Integ_MP2_BB,&
                              Lai_beta,Lai_beta,para_env)

      ! alpha-beta
      CALL contract_integrals(DI_AB,Emp2_RI_AB,DRI_AB,Emp2_AB*2.0_dp,dimen_RI,homo,homo_beta,virtual,virtual_beta,&
                              1.0_dp,1.0_dp,.FALSE.,&
                              Auto,Auto_beta,Integ_MP2_AB,&
                              Lai,Lai_beta,para_env)

      Emp2_RI=Emp2_RI_AA+Emp2_RI_BB+Emp2_RI_AB
      DRI=DRI_AA+DRI_BB+DRI_AB
      DI=DI_AA+DI_BB+DI_AB
    ELSE
      CALL contract_integrals(DI,Emp2_RI,DRI,Emp2,dimen_RI,homo,homo,virtual,virtual,&
                              2.0_dp,1.0_dp,.TRUE.,&
                              Auto,Auto,Integ_MP2,&
                              Lai,Lai,para_env)
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

  SUBROUTINE basis2p(nkind,ndof,basis_set_RI,RI_basis_parameter,p)
    INTEGER                                  :: nkind, ndof
    TYPE(gto_basis_set_type), DIMENSION(:), &
      POINTER                                :: basis_set_RI
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p

    INTEGER                                  :: ikind, ipos, iset

    p=0.0_dp
    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        p(ipos)=RI_basis_parameter(ikind)%zet(1,iset)
      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE p2basis(nkind,ndof,basis_set_RI,RI_basis_parameter,p)
    INTEGER                                  :: nkind, ndof
    TYPE(gto_basis_set_type), DIMENSION(:), &
      POINTER                                :: basis_set_RI
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: RI_basis_parameter
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p

    INTEGER                                  :: ikind, ipos, iset

    ipos=0
    DO ikind=1, nkind
      DO iset=1, RI_basis_parameter(ikind)%nset
        ipos=ipos+1
        RI_basis_parameter(ikind)%zet(1,iset)=p(ipos)
      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE contract_integrals(DI,Emp2_RI,DRI,Emp2,dimen_RI,homo,homo_beta,virtual,virtual_beta,&
                                fact,fact2,calc_ex,&
                                MOenerg,MOenerg_beta,abij,&
                                Lai,Lai_beta,para_env)
    REAL(KIND=dp)                            :: DI, Emp2_RI, DRI, Emp2
    INTEGER                                  :: dimen_RI, homo, homo_beta, &
                                                virtual, virtual_beta
    REAL(KIND=dp)                            :: fact, fact2
    LOGICAL                                  :: calc_ex
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: MOenerg, MOenerg_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: abij
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Lai, Lai_beta
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER                                  :: a, b, i, ij_counter, j
    REAL(KIND=dp)                            :: t_iajb, t_iajb_RI
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mat_ab

    ALLOCATE(mat_ab(virtual,virtual_beta))
  
    DI=0.0_dp
    Emp2_RI=0.0_dp
    ij_counter=0
    DO j=1, homo_beta
      DO i=1, homo
        ij_counter=ij_counter+1
        IF (MOD(ij_counter,para_env%num_pe)/=para_env%mepos) CYCLE
        mat_ab=0.0_dp
        mat_ab=MATMUL(TRANSPOSE(Lai(:,:,i)),Lai_beta(:,:,j))
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

  SUBROUTINE calc_elem_ij_proc(homo,homo_beta,para_env,elements_ij_proc,ij_list_proc)
    INTEGER                                  :: homo, homo_beta
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: elements_ij_proc
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_list_proc

    INTEGER                                  :: i, ij_counter, j

    elements_ij_proc=0
    ij_counter=-1
    DO i=1, homo
      DO j=1, homo_beta
        ij_counter=ij_counter+1
        IF (MOD(ij_counter,para_env%num_pe)==para_env%mepos) elements_ij_proc=elements_ij_proc+1
      END DO
    END DO

    ALLOCATE(ij_list_proc(elements_ij_proc,2))
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

END MODULE mp2_optimize_ri_basis

