!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_filter, cp_dbcsr_frobenius_norm, cp_dbcsr_get_occupation, &
       cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_nblkrows_total, &
       cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_external_control,             ONLY: external_control
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE dbcsr_types,                     ONLY: dbcsr_type_no_symmetry
  USE dm_ls_scf_methods,               ONLY: apply_matrix_preconditioner,&
                                             compute_matrix_preconditioner,&
                                             density_matrix_sign,&
                                             density_matrix_sign_fixed_mu
  USE dm_ls_scf_qs,                    ONLY: ls_scf_dm_to_ks,&
                                             ls_scf_init_qs,&
                                             matrix_ls_create,&
                                             matrix_ls_to_qs,&
                                             matrix_qs_to_ls
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
  USE f77_blas
  USE input_constants,                 ONLY: ls_cluster_atomic,&
                                             ls_cluster_molecular,&
                                             ls_s_inversion_hotelling,&
                                             ls_s_inversion_sign_sqrt,&
                                             ls_s_preconditioner_atomic,&
                                             ls_s_preconditioner_molecular,&
                                             ls_s_preconditioner_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                             matrix_sqrt_newton_schulz
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE molecule_types_new,              ONLY: molecule_of_atom,&
                                             molecule_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf'

  PUBLIC :: ls_scf

CONTAINS

! *****************************************************************************
!> \brief perform an linear scaling scf procedure: entry point 
!>       
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(ls_scf_env_type)                    :: ls_scf_env

    CALL timeset(routineN,handle)

    ! initialize the scf phase
    CALL ls_scf_init(qs_env,ls_scf_env,error)

    ! perform the actual scf
    CALL ls_scf_main(qs_env,ls_scf_env,error)

    ! do post scf processing
    CALL ls_scf_post(qs_env,ls_scf_env,error)

    CALL timestop(handle)

  END SUBROUTINE ls_scf

! *****************************************************************************
!> \brief Basic initialization of the LS type.
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! get basic quantities from the qs_env
    CALL get_qs_env(qs_env,nelectron_total=ls_scf_env%nelectron_total, &
                           matrix_s=matrix_s,&
                           dft_control=dft_control,&
                           molecule_set=molecule_set,&
                           input=input,&
                           has_unit_metric=ls_scf_env%has_unit_metric,&
                           para_env=ls_scf_env%para_env,&
                           nelectron_spin=ls_scf_env%nelectron_spin,error=error)

    ! copy some basic stuff
    ls_scf_env%nspins=dft_control%nspins
    ls_scf_env%natoms=cp_dbcsr_nblkrows_total(matrix_s(1)%matrix)
    CALL cp_para_env_retain(ls_scf_env%para_env,error)

    ! initialize block to group to defined molecules
    ALLOCATE(ls_scf_env%ls_mstruct%atom_to_molecule(ls_scf_env%natoms))
    CALL molecule_of_atom(molecule_set,atom_to_mol=ls_scf_env%ls_mstruct%atom_to_molecule,error=error)

    ! parse the ls_scf section and set derived quantities
    CALL ls_scf_init_read_write_input(input,ls_scf_env,unit_nr,error)

    ! some basic initialization of the QS side of things
    CALL ls_scf_init_qs(qs_env,error)

    ! create the matrix template for use in the ls procedures
    CALL matrix_ls_create(matrix_ls=ls_scf_env%matrix_s,matrix_qs=matrix_s(1)%matrix,& 
                          ls_mstruct=ls_scf_env%ls_mstruct,error=error)

    nspin=ls_scf_env%nspins
    ALLOCATE(ls_scf_env%matrix_p(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(ls_scf_env%matrix_p(ispin),error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_p(ispin),template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
    ENDDO

    ALLOCATE(ls_scf_env%matrix_ks(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(ls_scf_env%matrix_ks(ispin),error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_ks(ispin),template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
    ENDDO

    ! set up matrix S, and needed functions of S
    CALL ls_scf_init_matrix_s(matrix_s(1)%matrix,ls_scf_env,error)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init

! *****************************************************************************
!> \brief initialize S matrix related properties (sqrt, inverse...)
!>        Might be factored-out since this seems common code with the other SCF.
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init_matrix_S(matrix_s,ls_scf_env,error)
    TYPE(cp_dbcsr_type)                      :: matrix_s
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_matrix_S', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    REAL(KIND=dp)                            :: frob_matrix, frob_matrix_base
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1, matrix_tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! make our own copy of S
    IF (ls_scf_env%has_unit_metric) THEN
       CALL cp_dbcsr_set(ls_scf_env%matrix_s,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(ls_scf_env%matrix_s,1.0_dp,error=error)
    ELSE
       CALL matrix_qs_to_ls(ls_scf_env%matrix_s,matrix_s,ls_scf_env%ls_mstruct,error=error)
    ENDIF

    CALL cp_dbcsr_filter(ls_scf_env%matrix_s,ls_scf_env%eps_filter,error=error)

    ! needs a preconditioner for S
    IF (ls_scf_env%has_s_preconditioner) THEN
       CALL cp_dbcsr_init(ls_scf_env%matrix_bs_sqrt,error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_bs_sqrt,template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_init(ls_scf_env%matrix_bs_sqrt_inv,error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_bs_sqrt_inv,template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL compute_matrix_preconditioner(ls_scf_env%matrix_s,&
                            ls_scf_env%s_preconditioner_type, ls_scf_env%ls_mstruct, &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,&
                            ls_scf_env%eps_filter,error)
    ENDIF

    ! precondition S
    IF (ls_scf_env%has_s_preconditioner) THEN
       CALL apply_matrix_preconditioner(ls_scf_env%matrix_s,"forward", &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
    ENDIF

    ! compute sqrt(S) and inv(sqrt(S))
    IF (ls_scf_env%use_s_sqrt) THEN
        CALL cp_dbcsr_init(ls_scf_env%matrix_s_sqrt,error=error)
        CALL cp_dbcsr_init(ls_scf_env%matrix_s_sqrt_inv,error=error)
        CALL cp_dbcsr_create(ls_scf_env%matrix_s_sqrt,template=ls_scf_env%matrix_s,&
                             matrix_type=dbcsr_type_no_symmetry,error=error) 
        CALL cp_dbcsr_create(ls_scf_env%matrix_s_sqrt_inv,template=ls_scf_env%matrix_s,&
                             matrix_type=dbcsr_type_no_symmetry,error=error) 


        CALL matrix_sqrt_Newton_Schulz(ls_scf_env%matrix_s_sqrt,ls_scf_env%matrix_s_sqrt_inv,&
                                       ls_scf_env%matrix_s,ls_scf_env%eps_filter,error)

        
        IF (.TRUE.) THEN
           CALL cp_dbcsr_init(matrix_tmp1,error=error)
           CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,&
                                matrix_type=dbcsr_type_no_symmetry,error=error) 
           CALL cp_dbcsr_init(matrix_tmp2,error=error)
           CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,&
                                matrix_type=dbcsr_type_no_symmetry,error=error) 

           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_s,&
                                        0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)

           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                        0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)

           frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
           CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
           frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
           IF (unit_nr>0) THEN
              WRITE(unit_nr,*) "Error for (inv(sqrt(S))*S*inv(sqrt(S))-I)",frob_matrix/frob_matrix_base
           ENDIF

           CALL cp_dbcsr_release(matrix_tmp1,error=error) 
           CALL cp_dbcsr_release(matrix_tmp2,error=error) 
        ENDIF
    ENDIF

    ! compute the inverse of S
    IF (ls_scf_env%needs_s_inv) THEN
        CALL cp_dbcsr_init(ls_scf_env%matrix_s_inv,error=error)
        CALL cp_dbcsr_create(ls_scf_env%matrix_s_inv,template=ls_scf_env%matrix_s,&
                             matrix_type=dbcsr_type_no_symmetry,error=error) 
        IF (.NOT.ls_scf_env%use_s_sqrt) THEN
           CALL invert_Hotelling(ls_scf_env%matrix_s_inv,ls_scf_env%matrix_s,ls_scf_env%eps_filter,error)
        ELSE
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_s_sqrt_inv, &
                                        0.0_dp, ls_scf_env%matrix_s_inv, filter_eps=ls_scf_env%eps_filter,error=error)
        ENDIF
        IF (.TRUE.) THEN
           CALL cp_dbcsr_init(matrix_tmp1,error=error)
           CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,&
                                matrix_type=dbcsr_type_no_symmetry,error=error) 
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_inv, ls_scf_env%matrix_s,&
                                        0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
           frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
           CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
           frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
           IF (unit_nr>0) THEN
              WRITE(unit_nr,*) "Error for (inv(S)*S-I)",frob_matrix/frob_matrix_base
           ENDIF
           CALL cp_dbcsr_release(matrix_tmp1,error=error)
        ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_matrix_s

! *****************************************************************************
!> \brief parse the input section, no need to pass it around
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init_read_write_input(input,ls_scf_env,unit_nr,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(ls_scf_env_type), INTENT(INOUT)     :: ls_scf_env
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_read_write_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: mu
    TYPE(section_vals_type), POINTER         :: ls_scf_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ls_scf_section => section_vals_get_subs_vals(input,"DFT%LS_SCF",error=error)

    ! should come from input
    CALL section_vals_val_get(ls_scf_section,"EPS_SCF",r_val=ls_scf_env%eps_scf,error=error)
    CALL section_vals_val_get(ls_scf_section,"EPS_FILTER",r_val=ls_scf_env%eps_filter,error=error)
    CALL section_vals_val_get(ls_scf_section,"MU",r_val=mu,error=error)
    CALL section_vals_val_get(ls_scf_section,"FIXED_MU",l_val=ls_scf_env%fixed_mu,error=error)
    ls_scf_env%mu_spin=mu
    CALL section_vals_val_get(ls_scf_section,"MIXING_FRACTION",r_val=ls_scf_env%mixing_fraction,error=error)
    CALL section_vals_val_get(ls_scf_section,"MAX_SCF",i_val=ls_scf_env%max_scf,error=error)
    CALL section_vals_val_get(ls_scf_section,"S_PRECONDITIONER",i_val=ls_scf_env%s_preconditioner_type,error=error)
    CALL section_vals_val_get(ls_scf_section,"MATRIX_CLUSTER_TYPE",i_val=ls_scf_env%ls_mstruct%cluster_type,error=error)
    CALL section_vals_val_get(ls_scf_section,"SINGLE_PRECISION_MATRICES",l_val=ls_scf_env%ls_mstruct%single_precision,error=error)
    CALL section_vals_val_get(ls_scf_section,"S_INVERSION",i_val=ls_scf_env%s_inversion_type,error=error)
    CALL section_vals_val_get(ls_scf_section,"REPORT_ALL_SPARSITIES",l_val=ls_scf_env%report_all_sparsities,error=error)
    CALL section_vals_val_get(ls_scf_section,"PERFORM_MU_SCAN",l_val=ls_scf_env%perform_mu_scan,error=error)

    SELECT CASE(ls_scf_env%s_inversion_type) 
    CASE(ls_s_inversion_sign_sqrt)
        ls_scf_env%needs_s_inv=.TRUE.
        ls_scf_env%use_s_sqrt=.TRUE.
    CASE(ls_s_inversion_hotelling)
        ls_scf_env%needs_s_inv=.TRUE.
        ls_scf_env%use_s_sqrt=.FALSE.
    CASE DEFAULT
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    SELECT CASE(ls_scf_env%s_preconditioner_type)
    CASE(ls_s_preconditioner_none)
      ls_scf_env%has_s_preconditioner=.FALSE.
    CASE DEFAULT
      ls_scf_env%has_s_preconditioner=.TRUE.
    END SELECT

    ! an undocumented feature ... allows for just doing the initial guess, no expensive stuff
    IF (ls_scf_env%max_scf<0) THEN
        ls_scf_env%needs_s_inv=.FALSE.
        ls_scf_env%use_s_sqrt=.FALSE.
      ls_scf_env%has_s_preconditioner=.FALSE.
    ENDIF

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()')
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",30)," Linear scaling SCF ",REPEAT("-",29)
       WRITE(unit_nr,'(T2,A,T38,E20.3)') "eps_scf:",ls_scf_env%eps_scf
       WRITE(unit_nr,'(T2,A,T38,E20.3)') "eps_filter:",ls_scf_env%eps_filter
       WRITE(unit_nr,'(T2,A,T38,E20.3)') "mixing_fraction:",ls_scf_env%mixing_fraction
       WRITE(unit_nr,'(T2,A,T38,I20)') "max_scf:",ls_scf_env%max_scf
       WRITE(unit_nr,'(T2,A,T38,L20)') "fixed chemical potential (mu)",ls_scf_env%fixed_mu
       WRITE(unit_nr,'(T2,A,T38,L20)') "has unit metric:",ls_scf_env%has_unit_metric
       WRITE(unit_nr,'(T2,A,T38,L20)') "Computing inv(S):",ls_scf_env%needs_s_inv
       WRITE(unit_nr,'(T2,A,T38,L20)') "Computing sqrt(S):",ls_scf_env%use_s_sqrt
       WRITE(unit_nr,'(T2,A,T38,L20)') "Computing s preconditioner ",ls_scf_env%has_s_preconditioner
       WRITE(unit_nr,'(T2,A,T38,L20)') "Use single precision matrices",ls_scf_env%ls_mstruct%single_precision

       SELECT CASE(ls_scf_env%s_preconditioner_type)
       CASE(ls_s_preconditioner_none)
           WRITE(unit_nr,'(T2,A,T38,A20)') "S preconditioner type ","NONE"
       CASE(ls_s_preconditioner_atomic)
           WRITE(unit_nr,'(T2,A,T38,A20)') "S preconditioner type ","ATOMIC"
       CASE(ls_s_preconditioner_molecular)
           WRITE(unit_nr,'(T2,A,T38,A20)') "S preconditioner type ","MOLECULAR"
       END SELECT

       SELECT CASE(ls_scf_env%ls_mstruct%cluster_type)
       CASE(ls_cluster_atomic)
           WRITE(unit_nr,'(T2,A,T38,A20)') "Cluster type","ATOMIC"
       CASE(ls_cluster_molecular)
           WRITE(unit_nr,'(T2,A,T38,A20)') "Cluster type","MOLECULAR"
       END SELECT


       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
       WRITE(unit_nr,'()')
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_read_write_input

! *****************************************************************************
!> \brief Main SCF routine. Can we keep it clean ?
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_main(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_main', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iscf, ispin, &
                                                nelectron_spin_real, nspin, &
                                                unit_nr
    LOGICAL                                  :: should_stop
    REAL(KIND=dp)                            :: energy_diff, energy_new, &
                                                energy_old, t1, t2
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_mixing_old
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=ls_scf_env%nspins

    ! old quantities, useful for mixing
    ALLOCATE(matrix_mixing_old(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(matrix_mixing_old(ispin),error=error)
       CALL cp_dbcsr_create(matrix_mixing_old(ispin),template=ls_scf_env%matrix_ks(ispin), error=error)
    ENDDO

    energy_old=0.0_dp
    iscf=0

    ! the real SCF loop
    DO 

      ! check on max SCF or timing/exit
      CALL external_control(should_stop,"SCF",start_time=qs_env%start_time,target_time=qs_env%target_time,error=error)
      IF (should_stop .OR. iscf>=ls_scf_env%max_scf) THEN
         IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "SCF not converged! "
         EXIT
      ENDIF

      t1 = m_walltime()
      iscf=iscf+1

      ! first get a copy of the current KS matrix
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      DO ispin=1,nspin
         CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin),matrix_ks(ispin)%matrix,ls_scf_env%ls_mstruct,error=error)
         IF (ls_scf_env%has_s_preconditioner) THEN
             CALL apply_matrix_preconditioner(ls_scf_env%matrix_ks(ispin),"forward", &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
         ENDIF
         CALL cp_dbcsr_filter(ls_scf_env%matrix_ks(ispin),ls_scf_env%eps_filter,error=error)
      ENDDO

      ! initialize the mixing matrix with the current state if needed
      IF (iscf==1) THEN
         DO ispin=1,nspin
            CALL cp_dbcsr_copy(matrix_mixing_old(ispin), ls_scf_env%matrix_ks(ispin), error=error)
         ENDDO
      ENDIF

      ! perform the mixing of ks matrices
      DO ispin=1,nspin
         CALL cp_dbcsr_add(matrix_mixing_old(ispin), ls_scf_env%matrix_ks(ispin), &
                           1.0_dp-ls_scf_env%mixing_fraction,ls_scf_env%mixing_fraction, &
                           error=error)
      END DO

      ! compute the density matrix that matches it
      DO ispin=1,nspin
         ! we need the proper number of states
         nelectron_spin_real=ls_scf_env%nelectron_spin(ispin)
         IF (ls_scf_env%nspins==1) nelectron_spin_real=nelectron_spin_real/2

         CALL density_matrix_sign(ls_scf_env%matrix_p(ispin),ls_scf_env%mu_spin(ispin), ls_scf_env%fixed_mu, &
                                  matrix_mixing_old(ispin),ls_scf_env%matrix_s, &
                                  ls_scf_env%matrix_s_inv, &
                                  nelectron_spin_real,ls_scf_env%eps_filter,error)

         IF (ls_scf_env%has_s_preconditioner) THEN
             CALL apply_matrix_preconditioner(ls_scf_env%matrix_p(ispin),"forward", &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
         ENDIF
         CALL cp_dbcsr_filter(ls_scf_env%matrix_p(ispin),ls_scf_env%eps_filter,error=error)

         IF (ls_scf_env%nspins==1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin),2.0_dp,error=error)
      ENDDO

      ! compute the corresponding new energy KS matrix and new energy
      CALL ls_scf_dm_to_ks(qs_env,ls_scf_env,energy_new,error)

      ! report current SCF loop
      energy_diff=energy_new-energy_old
      energy_old=energy_new

      t2 = m_walltime()
      IF (unit_nr>0) THEN
         WRITE(unit_nr,*)
         WRITE(unit_nr,'(T2,A,I6,F20.9,F20.9,F12.6)') "SCF",iscf, energy_new,energy_diff, t2-t1
         WRITE(unit_nr,*)
      ENDIF

      ! exit criterium on the energy only for the time being
      IF (ABS(energy_diff)<ls_scf_env%eps_scf*ls_scf_env%nelectron_total) EXIT

    ENDDO

    ! free storage
    DO ispin=1,nspin
       CALL cp_dbcsr_release(matrix_mixing_old(ispin),error=error)
    ENDDO
    DEALLOCATE(matrix_mixing_old)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_main


! *****************************************************************************
!> \brief after SCF we have a density matrix, and the self consistent KS matrix
!>        analyze its properties.
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_post(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_post', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, unit_nr
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_w
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! compute the matrix_w if associated
    CALL get_qs_env(qs_env,matrix_w=matrix_w,error=error)
    IF (ASSOCIATED(matrix_w)) THEN
       CALL calculate_w_matrix(matrix_w,ls_scf_env,error=error)
    ENDIF

    ! compute properties

    IF (ls_scf_env%perform_mu_scan) CALL post_scf_mu_scan(ls_scf_env,error)

    IF (ls_scf_env%report_all_sparsities) CALL post_scf_sparsities(ls_scf_env,error)

    IF (.TRUE.) CALL post_scf_experiment(ls_scf_env,error)

    ! clean up used data

    CALL cp_dbcsr_release(ls_scf_env%matrix_s,error=error)

    IF (ls_scf_env%has_s_preconditioner) THEN
       CALL cp_dbcsr_release(ls_scf_env%matrix_bs_sqrt,error=error)
       CALL cp_dbcsr_release(ls_scf_env%matrix_bs_sqrt_inv,error=error)
    ENDIF

    IF (ls_scf_env%needs_s_inv) THEN
        CALL cp_dbcsr_release(ls_scf_env%matrix_s_inv,error=error)
    ENDIF

    IF (ls_scf_env%use_s_sqrt) THEN
        CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt,error=error)
        CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt_inv,error=error)
    ENDIF

    DO ispin=1,SIZE(ls_scf_env%matrix_p)
       CALL cp_dbcsr_release(ls_scf_env%matrix_p(ispin),error=error)
    ENDDO
    DEALLOCATE(ls_scf_env%matrix_p)

    DO ispin=1,SIZE(ls_scf_env%matrix_ks)
       CALL cp_dbcsr_release(ls_scf_env%matrix_ks(ispin),error=error)
    ENDDO
    DEALLOCATE(ls_scf_env%matrix_ks)

    CALL cp_para_env_release(ls_scf_env%para_env,error)

    DEALLOCATE(ls_scf_env%ls_mstruct%atom_to_molecule)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_post

! *****************************************************************************
!> \brief Compute the density matrix for various values of the chemical potential
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_mu_scan(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_mu_scan', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imu, ispin, &
                                                nelectron_spin_real, nmu, &
                                                nspin, unit_nr
    REAL(KIND=dp)                            :: mu, t1, t2, trace
    TYPE(cp_dbcsr_type)                      :: matrix_p
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=ls_scf_env%nspins
    
    CALL cp_dbcsr_init(matrix_p,error=error)
    CALL cp_dbcsr_create(matrix_p,template=ls_scf_env%matrix_p(1),error=error)

    nmu=10
    DO imu=0,nmu

       t1 = m_walltime()

       mu = -0.4_dp+imu*(0.1_dp+0.4_dp)/nmu

       IF (unit_nr>0) WRITE(unit_nr,*) "------- starting with mu ",mu

       DO ispin=1,nspin
          ! we need the proper number of states
          nelectron_spin_real=ls_scf_env%nelectron_spin(ispin)
          IF (ls_scf_env%nspins==1) nelectron_spin_real=nelectron_spin_real/2

          CALL density_matrix_sign_fixed_mu(matrix_p,trace,mu, &
                                            ls_scf_env%matrix_ks(ispin),ls_scf_env%matrix_s,&
                                            ls_scf_env%matrix_s_inv,ls_scf_env%eps_filter,error)
       ENDDO

       t2 = m_walltime()

       IF (unit_nr>0) WRITE(unit_nr,*) " obtained " ,mu,trace,t2-t1

    ENDDO

    CALL cp_dbcsr_release(matrix_p,error=error)

    CALL timestop(handle)

  END SUBROUTINE post_scf_mu_scan

! *****************************************************************************
!> \brief Report on the sparsities of various interesting matrices.
!>  
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_sparsities(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_sparsities', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: title
    INTEGER                                  :: handle, ispin, nspin, unit_nr
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1, matrix_tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=ls_scf_env%nspins

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()') 
       WRITE(unit_nr,'(T2,A,E17.3)') "Sparsity reports for eps_filter: ", ls_scf_env%eps_filter
       WRITE(unit_nr,'()') 
    ENDIF

    CALL report_matrix_sparsity(ls_scf_env%matrix_s,unit_nr,"overlap matrix (S)", &
                                ls_scf_env%eps_filter,error)

    DO ispin=1,nspin
       WRITE(title,'(A,I3)') "Kohn-Sham matrix (H) for spin ",ispin
       CALL report_matrix_sparsity(ls_scf_env%matrix_ks(ispin),unit_nr,title, &
                                   ls_scf_env%eps_filter,error)
    ENDDO

    CALL cp_dbcsr_init(matrix_tmp1,error=error)
    CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp2,error=error)
    CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)

    DO ispin=1,nspin
       WRITE(title,'(A,I3)') "Density matrix (P) for spin ",ispin
       CALL report_matrix_sparsity(ls_scf_env%matrix_p(ispin),unit_nr,title, &
                                   ls_scf_env%eps_filter,error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(ispin), &
                              0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)

       WRITE(title,'(A,I3,A)') "S * P(",ispin,")"
       CALL report_matrix_sparsity(matrix_tmp1,unit_nr,title, ls_scf_env%eps_filter,error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s, &
                              0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter, error=error)
       WRITE(title,'(A,I3,A)') "S * P(",ispin,") * S"
       CALL report_matrix_sparsity(matrix_tmp2,unit_nr,title, ls_scf_env%eps_filter,error)
    ENDDO

    IF (ls_scf_env%needs_s_inv) THEN
       CALL report_matrix_sparsity(ls_scf_env%matrix_s_inv,unit_nr,"inv(S)", &
                                   ls_scf_env%eps_filter,error)
       DO ispin=1,nspin
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_inv, ls_scf_env%matrix_ks(ispin), &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)

          WRITE(title,'(A,I3,A)') "inv(S) * H(",ispin,")"
          CALL report_matrix_sparsity(matrix_tmp1,unit_nr,title, ls_scf_env%eps_filter,error)
       ENDDO
    ENDIF

    IF (ls_scf_env%use_s_sqrt) THEN

       CALL report_matrix_sparsity(ls_scf_env%matrix_s_sqrt,unit_nr,"sqrt(S)", &
                                   ls_scf_env%eps_filter,error)
       CALL report_matrix_sparsity(ls_scf_env%matrix_s_sqrt_inv,unit_nr,"inv(sqrt(S))", &
                                   ls_scf_env%eps_filter,error)

       DO ispin=1,nspin
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin), &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                 0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)
          WRITE(title,'(A,I3,A)') "inv(sqrt(S)) * H(",ispin,") * inv(sqrt(S))"
          CALL report_matrix_sparsity(matrix_tmp2,unit_nr,title, ls_scf_env%eps_filter,error)
       ENDDO

       DO ispin=1,nspin
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%matrix_p(ispin), &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt, &
                                 0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)
          WRITE(title,'(A,I3,A)') "sqrt(S) * P(",ispin,") * sqrt(S)"
          CALL report_matrix_sparsity(matrix_tmp2,unit_nr,title, ls_scf_env%eps_filter,error)
       ENDDO

    ENDIF

    CALL cp_dbcsr_release(matrix_tmp1,error=error)
    CALL cp_dbcsr_release(matrix_tmp2,error=error)

    CALL timestop(handle)

  END SUBROUTINE post_scf_sparsities

! *****************************************************************************
!> \brief Helper routine to report on the sparsity of a single matrix,
!>        for several filtering values
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE report_matrix_sparsity(matrix,unit_nr,title,eps,error)
    TYPE(cp_dbcsr_type)                      :: matrix
    INTEGER                                  :: unit_nr
    CHARACTER(LEN=*)                         :: title
    REAL(KIND=dp)                            :: eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'report_matrix_sparsity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: eps_local, occ
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix,name=TRIM(title),error=error)
    CALL cp_dbcsr_copy(matrix_tmp,matrix,name=TRIM(title),error=error)

    IF (unit_nr>0) THEN
        WRITE(unit_nr,'(T2,A)') "Sparsity for : "//TRIM(title)
    ENDIF

    eps_local=eps
    DO
      IF (eps_local>1.1_dp) EXIT
      CALL cp_dbcsr_filter(matrix_tmp,eps_local,error=error)
      occ=cp_dbcsr_get_occupation(matrix_tmp)
      IF (unit_nr>0) WRITE(unit_nr,'(T2,F16.12,A3,F16.12)') eps_local," : ",occ
      eps_local=eps_local*10
    ENDDO

    CALL cp_dbcsr_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE report_matrix_sparsity

! *****************************************************************************
!> \brief Compute matrix_w as needed for the forces
!> \par History
!>       2010.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE calculate_w_matrix(matrix_w,ls_scf_env,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_w
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: scaling
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1, matrix_tmp2, &
                                                matrix_tmp3

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_init(matrix_tmp1,error=error)
    CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp2,error=error)
    CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp3,error=error)
    CALL cp_dbcsr_create(matrix_tmp3,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)

    IF (ls_scf_env%nspins==1) THEN
      scaling=0.5_dp
    ELSE
      scaling=1.0_dp
    ENDIF


    DO ispin=1,ls_scf_env%nspins

      CALL cp_dbcsr_copy(matrix_tmp3,ls_scf_env%matrix_ks(ispin),error=error)
      IF (ls_scf_env%has_s_preconditioner) THEN
         CALL apply_matrix_preconditioner(matrix_tmp3,"backward", &
                       ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
      ENDIF
      CALL cp_dbcsr_filter(matrix_tmp3,ls_scf_env%eps_filter,error=error)

      CALL cp_dbcsr_multiply("N", "N", scaling, ls_scf_env%matrix_p(ispin), matrix_tmp3, &
                              0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_p(ispin),&
                              0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)
      CALL matrix_ls_to_qs(matrix_w(ispin)%matrix, matrix_tmp2, ls_scf_env%ls_mstruct, error=error)
    ENDDO

    CALL cp_dbcsr_release(matrix_tmp1,error=error)
    CALL cp_dbcsr_release(matrix_tmp2,error=error)
    CALL cp_dbcsr_release(matrix_tmp3,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_w_matrix

! *****************************************************************************
!> \brief a place for quick experiments
!> \par History
!>       2010.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_experiment(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_experiment', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE post_scf_experiment

END MODULE dm_ls_scf
