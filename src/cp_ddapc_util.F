!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Density Derived atomic point charges from a QM calculation
!>      (see Bloechl, J. Chem. Phys. Vol. 103 pp. 7422-7428)
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp_ddapc_util
  USE atomic_charges,                  ONLY: print_atomic_charges
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: ddapc_restraint_type,&
                                             dft_control_type
  USE cp_ddapc_forces,                 ONLY: evaluate_restraint_functional
  USE cp_ddapc_methods,                ONLY: build_A_matrix,&
                                             build_b_vector,&
                                             build_der_A_matrix_rows,&
                                             build_der_b_vector,&
                                             cleanup_g_dot_rvec_sin_cos,&
                                             prep_g_dot_rvec_sin_cos
  USE cp_ddapc_types,                  ONLY: cp_ddapc_create,&
                                             cp_ddapc_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: do_full_density,&
                                             do_spin_density
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_transfer
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_util'
  PUBLIC :: get_ddapc,&
            restraint_functional_potential,&
            modify_hartree_pot,&
            cp_ddapc_init

CONTAINS

! *****************************************************************************
!> \brief Initialize the cp_ddapc_environment
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE cp_ddapc_init(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, iw2, &
                                                n_rep_val, num_gauss, stat
    LOGICAL                                  :: allocate_ddapc_env, failure, &
                                                unimplemented
    REAL(KIND=dp)                            :: gcut, pfact, rcmin, Vol
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_radii, radii
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_g
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: density_fit_section

    CALL timeset(routineN,handle)
    failure         = .FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(dft_control, rho, rho_tot_g, rho_core, rho0_s_gs, pw_env,&
            radii, inp_radii, particle_set, qs_charges)
    allocate_ddapc_env = qs_env%cp_ddapc_ewald%do_solvation.OR.&
                         qs_env%cp_ddapc_ewald%do_qmmm_periodic_decpl.OR.&
                         qs_env%cp_ddapc_ewald%do_decoupling.OR.&
                         qs_env%cp_ddapc_ewald%do_restraint
    unimplemented =      qs_env%dft_control%qs_control%semi_empirical.OR.&
                         qs_env%dft_control%qs_control%dftb
    IF (allocate_ddapc_env.AND.unimplemented) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="DDAP charges work only with GPW/GAPW code.", &
            error=error, error_level=cp_failure_level)
    END IF
    allocate_ddapc_env = allocate_ddapc_env .OR.&
                         qs_env%cp_ddapc_ewald%do_property
    allocate_ddapc_env = allocate_ddapc_env .AND.(.NOT.unimplemented)
    IF (.NOT.failure.AND.allocate_ddapc_env) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       dft_control=dft_control,&
                       rho=rho,&
                       rho_core=rho_core,&
                       rho0_s_gs=rho0_s_gs,&
                       pw_env=pw_env,&
                       qs_charges=qs_charges,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)
       density_fit_section => section_vals_get_subs_vals(qs_env%input,"DFT%DENSITY_FITTING", error=error)
       iw=cp_print_key_unit_nr(logger,density_fit_section,&
            "PROGRAM_RUN_INFO",".FitCharge",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(/,A)')" Initializing the DDAPC Environment"
       END IF
       CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pool,error=error)
       CALL pw_pool_create_pw(auxbas_pool,rho_tot_g, in_space=RECIPROCALSPACE,&
            use_data=COMPLEXDATA1D, error=error)
       Vol = rho_tot_g%pw_grid%vol
       !
       ! Get Input Parameters
       !
       CALL section_vals_val_get(density_fit_section,"RADII",n_rep_val=n_rep_val,error=error)
       IF (n_rep_val /= 0) THEN
          CALL section_vals_val_get(density_fit_section,"RADII",r_vals=inp_radii,error=error)
          num_gauss = SIZE(inp_radii)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          radii = inp_radii
       ELSE
          CALL section_vals_val_get(density_fit_section,"NUM_GAUSS",i_val=num_gauss,error=error)
          CALL section_vals_val_get(density_fit_section,"MIN_RADIUS",r_val=rcmin,error=error)
          CALL section_vals_val_get(density_fit_section,"PFACTOR",r_val=pfact,error=error)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, num_gauss
             radii(i) = rcmin*pfact**(i-1)
          END DO
       END IF
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       ! Create DDAPC environment
       iw2=cp_print_key_unit_nr(logger,density_fit_section,&
            "PROGRAM_RUN_INFO/CONDITION_NUMBER",".FitCharge",error=error)
       ! Initialization of the cp_ddapc_env and of the cp_ddapc_ewald environment
       !NB pass qs_env%para_env for parallelization of ewald_ddapc_pot()
       CALL cp_ddapc_create(qs_env%para_env,&
                            qs_env%cp_ddapc_env,&
                            qs_env%cp_ddapc_ewald,&
                            particle_set,&
                            radii,&
                            cell,&
                            super_cell,&
                            rho_tot_g,&
                            gcut,&
                            iw2,&
                            Vol, &
                            qs_env%input,&
                            error)
       CALL cp_print_key_finished_output(iw2,logger,density_fit_section,&
            "PROGRAM_RUN_INFO/CONDITION_NUMBER",error=error)
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_g,&
            accept_non_compatible=.TRUE.,error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_ddapc_init

! *****************************************************************************
!> \brief Computes the Density Derived Atomic Point Charges
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE get_ddapc(qs_env, calc_force, density_fit_section,&
       density_type,qout1, qout2, out_radii, dq_out, ext_rho_tot_g,&
       Itype_of_density,iwc, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    TYPE(section_vals_type), POINTER         :: density_fit_section
    INTEGER, OPTIONAL                        :: density_type
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: qout1, qout2, out_radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dq_out
    TYPE(pw_type), OPTIONAL, POINTER         :: ext_rho_tot_g
    CHARACTER(LEN=*), OPTIONAL               :: Itype_of_density
    INTEGER, INTENT(IN), OPTIONAL            :: iwc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_ddapc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: type_of_density
    INTEGER :: handle, handle2, handle3, i, ii, iparticle, iparticle0, ispin, &
      iw, j, myid, n_rep_val, ndim, nparticles, num_gauss, pmax, pmin, stat
    LOGICAL                                  :: failure, need_f
    REAL(KIND=dp)                            :: c1, c3, ch_dens, gcut, pfact, &
                                                rcmin, Vol
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: AmI_bv, AmI_cv, bv, cv, &
                                                cvT_AmI, cvT_AmI_dAmj, &
                                                dAmj_qv, qtot, qv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dbv, g_dot_rvec_cos, &
                                                g_dot_rvec_sin
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: dAm, dqv, tv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_radii, radii
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_g
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_rho_type), POINTER               :: rho

!NB variables for doing build_der_A_matrix_rows in blocks
!NB refactor math in inner loop - no need for dqv0
!!NB refactor math in inner loop - new temporaries

    EXTERNAL dgemv, dgemm

    CALL timeset(routineN,handle)
    failure         = .FALSE.
    need_f          = .FALSE.
    IF (PRESENT(calc_force))  need_f      = calc_force
    logger => cp_error_get_logger(error)
    NULLIFY(dft_control, rho, rho_tot_g, rho_core, rho0_s_gs, pw_env,&
            radii, inp_radii, particle_set, qs_charges, cp_ddapc_env)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       dft_control=dft_control,&
                       rho=rho,&
                       rho_core=rho_core,&
                       rho0_s_gs=rho0_s_gs,&
                       pw_env=pw_env,&
                       qs_charges=qs_charges,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)
       IF (PRESENT(iwc)) THEN
          iw = iwc
       ELSE
          iw=cp_print_key_unit_nr(logger,density_fit_section,&
               "PROGRAM_RUN_INFO",".FitCharge",error=error)
       END IF
       CALL pw_env_get(pw_env=pw_env,&
                       auxbas_pw_pool=auxbas_pool,&
                       error=error)
       CALL pw_pool_create_pw(auxbas_pool,rho_tot_g, in_space=RECIPROCALSPACE,&
            use_data=COMPLEXDATA1D, error=error)
       IF (PRESENT(ext_rho_tot_g)) THEN
          ! If provided use the input density in g-space
          CALL pw_transfer(ext_rho_tot_g,rho_tot_g, error=error)
          type_of_density=Itype_of_density
       ELSE
          IF(PRESENT(density_type))THEN
             myid=density_type
          ELSE
             CALL section_vals_val_get(qs_env%input,&
                                    "PROPERTIES%FIT_CHARGE%TYPE_OF_DENSITY",i_val=myid,error=error)
          END IF
          SELECT CASE(myid)
          CASE(do_full_density)
             ! Otherwise build the total QS density (electron+nuclei) in G-space
             IF(dft_control%qs_control%gapw) THEN
                CALL pw_transfer(rho0_s_gs%pw,rho_tot_g, error=error)
             ELSE
                CALL pw_transfer(rho_core%pw,rho_tot_g, error=error)
             END IF
             DO ispin=1,SIZE(rho%rho_g)
                CALL pw_axpy(rho%rho_g(ispin)%pw,rho_tot_g, error=error)
             END DO
             type_of_density = "FULL DENSITY"
          CASE(do_spin_density)
             CALL pw_copy (rho%rho_g(1)%pw,rho_tot_g, error=error)
             CALL pw_axpy (rho%rho_g(2)%pw,rho_tot_g,alpha=-1._dp, error=error)
             type_of_density = "SPIN DENSITY"
          END SELECT
       END IF
       Vol = rho%rho_r(1)%pw%pw_grid%vol
       ch_dens = 0.0_dp
       IF (rho_tot_g%pw_grid%have_g0) ch_dens = rho_tot_g%cc(1)
       CALL mp_sum(ch_dens,logger%para_env%group)
       !
       ! Get Input Parameters
       !
       CALL section_vals_val_get(density_fit_section,"RADII",n_rep_val=n_rep_val,error=error)
       IF (n_rep_val /= 0) THEN
          CALL section_vals_val_get(density_fit_section,"RADII",r_vals=inp_radii,error=error)
          num_gauss = SIZE(inp_radii)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          radii = inp_radii
       ELSE
          CALL section_vals_val_get(density_fit_section,"NUM_GAUSS",i_val=num_gauss,error=error)
          CALL section_vals_val_get(density_fit_section,"MIN_RADIUS",r_val=rcmin,error=error)
          CALL section_vals_val_get(density_fit_section,"PFACTOR",r_val=pfact,error=error)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, num_gauss
             radii(i) = rcmin*pfact**(i-1)
          END DO
       END IF
       IF (PRESENT(out_radii)) THEN
          IF (ASSOCIATED(out_radii)) THEN
             DEALLOCATE(out_radii, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ALLOCATE(out_radii(SIZE(radii)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          out_radii = radii
       END IF
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       cp_ddapc_env => qs_env%cp_ddapc_env
       !
       ! Start with the linear system
       !
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(bv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(qv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(qtot(SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(cv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL timeset(routineN//"-charges",handle2)
       bv  = 0.0_dp
       cv  = 1.0_dp/Vol
       CALL build_b_vector(bv, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
            particle_set, radii, rho_tot_g, gcut, error);  bv = bv / Vol
       CALL mp_sum(bv,rho_tot_g%pw_grid%para%group)
       c1 = DOT_PRODUCT(cv,MATMUL(cp_ddapc_env%AmI,bv))-ch_dens
       c1 = c1 / cp_ddapc_env%c0
       qv = - MATMUL(cp_ddapc_env%AmI,(bv-c1*cv))
       j  = 0
       qtot = 0.0_dp
       DO i = 1, ndim, num_gauss
          j = j+1
          DO ii = 1, num_gauss
             qtot(j) = qtot(j) + qv((i-1)+ii)
          END DO
       END DO
       IF (PRESENT(qout1)) THEN
          IF (ASSOCIATED(qout1)) THEN
             CPPostcondition(SIZE(qout1)==SIZE(qv),cp_failure_level,routineP,error,failure)
          ELSE
             ALLOCATE(qout1(SIZE(qv)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          qout1 = qv
       END IF
       IF (PRESENT(qout2)) THEN
          IF (ASSOCIATED(qout2)) THEN
             CPPostcondition(SIZE(qout2)==SIZE(qtot),cp_failure_level,routineP,error,failure)
          ELSE
             ALLOCATE(qout2(SIZE(qtot)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          qout2 = qtot
       END IF
       CALL print_atomic_charges(particle_set,iw,title=" DDAP "//&
            TRIM(type_of_density)//" charges:",atomic_charges=qtot)
       CALL timestop(handle2)
       !
       ! If requested evaluate also the correction to derivatives due to Pulay Forces
       !
       IF (need_f) THEN
          CALL timeset(routineN//"-forces",handle3)
          IF (iw>0) THEN
             WRITE(iw,'(T3,A)')" Evaluating DDAPC atomic derivatives .."
          END IF
          ALLOCATE(dAm(ndim, ndim,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dbv(ndim,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dqv(ndim,SIZE(particle_set),3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !NB refactor math in inner loop - no more dqv0, but new temporaries instead
          ALLOCATE(cvT_AmI(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(cvT_AmI_dAmj(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(tv(ndim,SIZE(particle_set),3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(AmI_cv(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          cvT_AmI = MATMUL(cv, cp_ddapc_env%AmI)
          AmI_cv = MATMUL(cp_ddapc_env%AmI,cv)
          ALLOCATE(dAmj_qv(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(AmI_bv(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          AmI_bv = MATMUL(cp_ddapc_env%AmI,bv)

          !NB call routine to precompute sin(g.r) and cos(g.r),
          ! so it doesn't have to be done for each r_i-r_j pair in build_der_A_matrix_rows()
          CALL prep_g_dot_rvec_sin_cos(rho_tot_g, particle_set, gcut, g_dot_rvec_sin, g_dot_rvec_cos, error)
          !NB do build_der_A_matrix_rows in blocks, for more efficient use of DGEMM
#define NPSET 100
          DO iparticle0 = 1, SIZE(particle_set), NPSET
             nparticles = MIN(NPSET, SIZE(particle_set)-iparticle0+1)
             !NB each dAm is supposed to have one block of rows and one block of columns
             !NB for derivatives with respect to each atom.  build_der_A_matrix_rows()
             !NB just returns rows, since dAm is symmetric, and missing columns can be 
             !NB reconstructed with a simple transpose, as below
             CALL build_der_A_matrix_rows(dAm, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
                  particle_set, radii, rho_tot_g, gcut, iparticle0, nparticles, g_dot_rvec_sin, g_dot_rvec_cos, error); 
             !NB no more reduction of dbv and dAm - instead we go through with each node's contribution
             !NB and reduce resulting charges/forces once, at the end.  Intermediate speedup can be
             !NB had by reducing dqv after the inner loop, and then other routines don't need to know
             !NB that contributions to dqv are distributed over the nodes.
             !NB also get rid of zeroing of dAm and division by Vol**2 - it's slow, and can be done
             !NB more quickly later, to a scalar or vector rather than a matrix
             DO iparticle = iparticle0, iparticle0+nparticles-1
             IF (debug_this_module) THEN
                CALL debug_der_A_matrix(dAm, particle_set, radii, rho_tot_g,&
                     gcut, iparticle, Vol, qs_env, error)
                cp_ddapc_env => qs_env%cp_ddapc_env
             END IF
                dbv = 0.0_dp
                CALL build_der_b_vector(dbv, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
                     particle_set, radii, rho_tot_g, gcut, iparticle, error);  dbv = dbv / Vol
                IF (debug_this_module) THEN
                   CALL debug_der_b_vector(dbv, particle_set, radii, rho_tot_g,&
                        gcut, iparticle, Vol, qs_env, error)
                   cp_ddapc_env => qs_env%cp_ddapc_env 
                END IF
             DO j = 1, 3
                !NB dAmj is actually pretty sparse - one block of cols + one block of rows - use this here:
                pmin = (iparticle-1)*SIZE(radii)+1
                pmax = iparticle*SIZE(radii)
                   !NB multiply by block of columns that aren't explicitly in dAm, but can be reconstructured
                   !NB as transpose of relevant block of rows
                IF (pmin > 1) THEN
                     dAmj_qv(:pmin-1)      = MATMUL(TRANSPOSE(dAm(pmin:pmax,:pmin-1,j)),qv(pmin:pmax))
                     cvT_AmI_dAmj(:pmin-1) = MATMUL(TRANSPOSE(dAm(pmin:pmax,:pmin-1,j)),cvT_AmI(pmin:pmax))
                ENDIF
                   !NB multiply by block of rows that are explicitly in dAm
                dAmj_qv(pmin:pmax)      = MATMUL(dAm(pmin:pmax,:,j),qv(:))
                cvT_AmI_dAmj(pmin:pmax) = MATMUL(dAm(pmin:pmax,:,j),cvT_AmI(:))
                   !NB multiply by block of columns that aren't explicitly in dAm, but can be reconstructured
                   !NB as transpose of relevant block of rows
                IF (pmax < SIZE(particle_set)*SIZE(radii)) THEN
                     dAmj_qv(pmax+1:)      = MATMUL(TRANSPOSE(dAm(pmin:pmax,pmax+1:,j)),qv(pmin:pmax))
                     cvT_AmI_dAmj(pmax+1:) = MATMUL(TRANSPOSE(dAm(pmin:pmax,pmax+1:,j)),cvT_AmI(pmin:pmax))
                ENDIF
                dAmj_qv      = dAmj_qv / (Vol*Vol)
                cvT_AmI_dAmj = cvT_AmI_dAmj / (Vol*Vol)
                c3 = DOT_PRODUCT(cvT_AmI_dAmj,AmI_bv) - DOT_PRODUCT(cvT_AmI, dbv(:,j)) - c1*DOT_PRODUCT(cvT_AmI_dAmj,AmI_cv)
                tv(:,iparticle,j) = -(dAmj_qv + dbv(:,j) + c3/cp_ddapc_env%c0*cv)
             END DO ! j
             !NB zero relevant parts of dAm here
             dAm((iparticle-1)*SIZE(radii)+1:iparticle*SIZE(radii),:,:) = 0.0_dp
                !! dAm(:,(iparticle-1)*SIZE(radii)+1:iparticle*SIZE(radii),:) = 0.0_dp
          END DO ! iparticle
          END DO ! iparticle0
          !NB final part of refactoring of math - one dgemm is faster than many dgemv
          CALL dgemm('N','N',SIZE(dqv,1),SIZE(dqv,2)*SIZE(dqv,3),SIZE(cp_ddapc_env%AmI,2),1.0_dp,&
                     cp_ddapc_env%AmI,SIZE(cp_ddapc_env%AmI,1),tv,SIZE(tv,1),0.0_dp,dqv,SIZE(dqv,1))
          !NB deallocate g_dot_rvec_sin and g_dot_rvec_cos
          CALL cleanup_g_dot_rvec_sin_cos(g_dot_rvec_sin, g_dot_rvec_cos, error)
          !NB moved reduction out to where dqv is used to compute
          !NB  a force contribution (smaller array to reduce, just size(particle_set) x 3)
          !NB namely ewald_ddapc_force(), cp_decl_ddapc_forces(), restraint_functional_force()
          CPPostcondition(PRESENT(dq_out),cp_failure_level,routineP,error,failure)
          IF (.NOT.ASSOCIATED(dq_out)) THEN
             ALLOCATE(dq_out(SIZE(dqv,1),SIZE(dqv,2),SIZE(dqv,3)), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ELSE
             CPPostcondition(SIZE(dqv,1)==SIZE(dq_out,1),cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(dqv,2)==SIZE(dq_out,2),cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(dqv,3)==SIZE(dq_out,3),cp_failure_level,routineP,error,failure)
          END IF
          dq_out = dqv
          IF (debug_this_module) THEN
             CALL debug_charge(dqv, qs_env, density_fit_section, &
                  particle_set, radii, rho_tot_g, type_of_density, error)
             cp_ddapc_env => qs_env%cp_ddapc_env
          END IF
          DEALLOCATE(dqv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dAm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dbv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !NB deallocate new temporaries
          DEALLOCATE(cvT_AmI, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(cvT_AmI_dAmj, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(AmI_cv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(tv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dAmj_qv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(AmI_bv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL timestop(handle3)
       END IF
       !
       ! End of charge fit
       !
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(bv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qtot, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.PRESENT(iwc)) THEN
          CALL cp_print_key_finished_output(iw,logger,density_fit_section,&
               "PROGRAM_RUN_INFO",error=error)
       END IF
       CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_g,&
            accept_non_compatible=.TRUE.,error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE get_ddapc

! *****************************************************************************
!> \brief modify hartree potential to handle restraints in DDAPC scheme
!> \par History
!>      02.2006  modified [Teo]
! *****************************************************************************
SUBROUTINE restraint_functional_potential(v_hartree_gspace,&
     density_fit_section, particle_set, AmI, radii, charges,&
     ddapc_restraint_control,energy_res, error)
    TYPE(pw_p_type)                          :: v_hartree_gspace
    TYPE(section_vals_type), POINTER         :: density_fit_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: AmI
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii, charges
    TYPE(ddapc_restraint_type), &
      INTENT(INOUT)                          :: ddapc_restraint_control
    REAL(KIND=dp), INTENT(INOUT)             :: energy_res
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'restraint_functional_potential', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: g_corr, phase
    INTEGER                                  :: handle, idim, ig, igauss, &
                                                iparticle, n_gauss, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, fac, fac2, g2, gcut, &
                                                gcut2, gfunc, gvec(3), rc, &
                                                rc2, rvec(3), sfac, Vol, w
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cv, uv
    TYPE(pw_type), POINTER                   :: g_hartree

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(g_hartree)
    IF (.NOT.failure) THEN
       n_gauss = SIZE(radii)
       ALLOCATE(cv(n_gauss*SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(uv(n_gauss*SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       uv      = 0.0_dp
       CALL evaluate_restraint_functional(ddapc_restraint_control, n_gauss, uv,&
            charges, energy_res, error)
       !
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       gcut2 = gcut*gcut
       g_hartree =>  v_hartree_gspace%pw
       Vol = g_hartree%pw_grid%vol
       cv   =  1.0_dp/Vol
       sfac = -1.0_dp/Vol
       fac  = DOT_PRODUCT(cv, MATMUL(AmI,cv))
       fac2 = DOT_PRODUCT(cv, MATMUL(AmI,uv))
       cv   = uv - cv * fac2/fac
       cv   = MATMUL(AmI,cv)
       IF (g_hartree%pw_grid%have_g0) g_hartree%cc(1) = g_hartree%cc(1) + sfac*fac2/fac
       DO ig = g_hartree % pw_grid % first_gne0, g_hartree % pw_grid % ngpts_cut_local
          g2 = g_hartree % pw_grid % gsq ( ig )
          w  = 4.0_dp * pi * (g2 - gcut2)**2.0_dp/(g2*gcut2)
          IF (g2 > gcut2) EXIT
          gvec  = g_hartree % pw_grid % g (:,ig)
          g_corr = 0.0_dp
          idim = 0
          DO iparticle = 1, SIZE(particle_set)
             DO igauss = 1, SIZE(radii)
                idim   = idim + 1
                rc     = radii(igauss)
                rc2    = rc*rc
                rvec   = particle_set(iparticle)%r
                arg    = DOT_PRODUCT(gvec,rvec)
                phase  = CMPLX(COS(arg),-SIN(arg),KIND=dp)
                gfunc  = EXP(-g2*rc2/4.0_dp)
                g_corr = g_corr + gfunc * cv(idim) * phase
             END DO
          END DO
          g_corr = g_corr * w
          g_hartree%cc(ig) = g_hartree%cc(ig) + sfac*g_corr/Vol
       END DO
    END IF
    CALL timestop(handle)
  END SUBROUTINE restraint_functional_potential

! *****************************************************************************
!> \brief Modify the Hartree potential
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE modify_hartree_pot(v_hartree_gspace, density_fit_section,&
       particle_set, M, AmI, radii, charges, error)
    TYPE(pw_p_type)                          :: v_hartree_gspace
    TYPE(section_vals_type), POINTER         :: density_fit_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M, AmI
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii, charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'modify_hartree_pot', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: g_corr, phase
    INTEGER                                  :: handle, idim, ig, igauss, &
                                                iparticle, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: arg, fac, fac2, g2, gcut, &
                                                gcut2, gfunc, gvec(3), rc, &
                                                rc2, rvec(3), sfac, Vol, w
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cv, uv
    TYPE(pw_type), POINTER                   :: g_hartree

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(g_hartree)
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       gcut2 = gcut*gcut
       g_hartree =>  v_hartree_gspace%pw
       Vol = g_hartree%pw_grid%vol
       ALLOCATE(cv(SIZE(M,1)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(uv(SIZE(M,1)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       cv   = 1.0_dp/Vol
       uv   = MATMUL(M,charges)
       sfac = -1.0_dp/Vol
       fac  = DOT_PRODUCT(cv, MATMUL(AmI,cv))
       fac2 = DOT_PRODUCT(cv, MATMUL(AmI,uv))
       cv   = uv - cv * fac2/fac
       cv   = MATMUL(AmI,cv)
       IF (g_hartree%pw_grid%have_g0) g_hartree%cc(1) = g_hartree%cc(1) + sfac*fac2/fac
       DO ig = g_hartree % pw_grid % first_gne0, g_hartree % pw_grid % ngpts_cut_local
          g2 = g_hartree % pw_grid % gsq ( ig )
          w  = 4.0_dp * pi * (g2 - gcut2)**2.0_dp/(g2*gcut2)
          IF (g2 > gcut2) EXIT
          gvec  = g_hartree % pw_grid % g (:,ig)
          g_corr = 0.0_dp
          idim = 0
          DO iparticle = 1, SIZE(particle_set)
             DO igauss = 1, SIZE(radii)
                idim   = idim + 1
                rc     = radii(igauss)
                rc2    = rc*rc
                rvec   = particle_set(iparticle)%r
                arg    = DOT_PRODUCT(gvec,rvec)
                phase  = CMPLX(COS(arg),-SIN(arg),KIND=dp)
                gfunc  = EXP(-g2*rc2/4.0_dp)
                g_corr = g_corr + gfunc * cv(idim) * phase
             END DO
          END DO
          g_corr = g_corr * w
          g_hartree%cc(ig) = g_hartree%cc(ig) + sfac*g_corr/Vol
       END DO
    END IF
    CALL timestop(handle)
  END SUBROUTINE modify_hartree_pot

! *****************************************************************************
!> \brief To Debug the derivative of the B vector for the solution of the
!>      linear system
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE debug_der_b_vector(dbv,  particle_set, radii,&
       rho_tot_g, gcut, iparticle, Vol, qs_env, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: dbv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(in)                      :: iparticle
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_der_b_vector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, kk, ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, rvec(3), v0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: bv1, bv2, ddbv
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env

    NULLIFY(cp_ddapc_env)
    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       dx   = 0.01_dp
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(bv1(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(bv2(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ddbv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       rvec = particle_set(iparticle)%r
       cp_ddapc_env => qs_env%cp_ddapc_env
       DO i = 1, 3
          bv1 = 0.0_dp
          bv2 = 0.0_dp
          particle_set(iparticle)%r(i) = rvec(i) + dx
          CALL build_b_vector(bv1, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, error);  bv1 = bv1 / Vol
          CALL mp_sum(bv1,rho_tot_g%pw_grid%para%group)
          particle_set(iparticle)%r(i) = rvec(i) - dx
          CALL build_b_vector(bv2, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, error);  bv2 = bv2 / Vol
          CALL mp_sum(bv2,rho_tot_g%pw_grid%para%group)
          ddbv = (bv1 - bv2)/(2.0_dp*dx)
          DO kk = 1, SIZE(ddbv)
             IF (ddbv(kk).GT.1.0E-8_dp) THEN
                v0 = ABS(dbv(kk,i)-ddbv(kk))/ddbv(kk)*100.0_dp
                WRITE(*,*)"Error % on B ::",v0
                IF (v0.GT.0.1_dp) THEN
                   WRITE(*,'(A,2I5,2F15.9)')"ERROR IN DERIVATIVE OF B VECTOR, IPARTICLE, ICOORD:",iparticle,i,&
                        dbv(kk,i),ddbv(kk)
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END IF
          END DO
          particle_set(iparticle)%r = rvec
       END DO
       DEALLOCATE(bv1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(bv2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ddbv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE debug_der_b_vector

! *****************************************************************************
!> \brief To Debug the derivative of the A matrix for the solution of the
!>      linear system
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE debug_der_A_matrix(dAm, particle_set, radii,&
       rho_tot_g, gcut, iparticle, Vol, qs_env, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: dAm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(in)                      :: iparticle
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_der_A_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, kk, ll, ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, rvec(3), v0
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Am1, Am2, ddAm, &
                                                g_dot_rvec_cos, g_dot_rvec_sin
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env

!NB new temporaries sin(g.r) and cos(g.r), as used in get_ddapc, to speed up build_der_A_matrix()

    failure = .FALSE.
    NULLIFY(cp_ddapc_env)
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       dx   = 0.01_dp
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(Am1(ndim,ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Am2(ndim,ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ddAm(ndim,ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       rvec = particle_set(iparticle)%r
       cp_ddapc_env => qs_env%cp_ddapc_env
       CALL prep_g_dot_rvec_sin_cos(rho_tot_g, particle_set, gcut, g_dot_rvec_sin, g_dot_rvec_cos, error)
       DO i = 1, 3
          Am1 = 0.0_dp
          Am2 = 0.0_dp
          particle_set(iparticle)%r(i) = rvec(i) + dx
          CALL build_A_matrix(Am1, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, g_dot_rvec_sin, g_dot_rvec_cos, &
               error);  Am1 = Am1 /(Vol*Vol)
          CALL mp_sum(Am1,rho_tot_g%pw_grid%para%group)
          particle_set(iparticle)%r(i) = rvec(i) - dx
          CALL build_A_matrix(Am2, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, g_dot_rvec_sin, g_dot_rvec_cos, &
               error);  Am2 = Am2 /(Vol*Vol)
          CALL mp_sum(Am2,rho_tot_g%pw_grid%para%group)
          ddAm = (Am1 - Am2)/(2.0_dp*dx)
          DO kk = 1, SIZE(ddAm,1)
             DO ll = 1, SIZE(ddAm,2)
                IF (ddAm(kk,ll).GT.1.0E-8_dp) THEN
                   v0 = ABS(dAm(kk,ll,i)-ddAm(kk,ll))/ddAm(kk,ll)*100.0_dp
                   WRITE(*,*)"Error % on A ::",v0,Am1(kk,ll), Am2(kk,ll),iparticle,i,kk,ll
                   IF (v0.GT.0.1_dp) THEN
                      WRITE(*,'(A,4I5,2F15.9)')"ERROR IN DERIVATIVE OF A MATRIX, IPARTICLE, ICOORD:",iparticle,i,kk,ll,&
                           dAm(kk,ll,i),ddAm(kk,ll)
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
          END DO
          particle_set(iparticle)%r = rvec
       END DO
       CALL cleanup_g_dot_rvec_sin_cos(g_dot_rvec_sin, g_dot_rvec_cos, error)
       DEALLOCATE(Am1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(Am2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ddAm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE debug_der_A_matrix

! *****************************************************************************
!> \brief To Debug the fitted charges
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE debug_charge(dqv, qs_env, density_fit_section,&
       particle_set, radii, rho_tot_g, type_of_density, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: dqv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: density_fit_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    CHARACTER(LEN=*)                         :: type_of_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_charge', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iparticle, kk, &
                                                ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, rvec(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ddqv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: qtot1, qtot2

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       WRITE(*,*)"DEBUG_CHARGE_ROUTINE"
       ndim = SIZE(particle_set)*SIZE(radii)
       NULLIFY(qtot1, qtot2)
       ALLOCATE(qtot1(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(qtot2(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ddqv(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !
       dx = 0.001_dp
       DO iparticle = 1, SIZE(particle_set)
          rvec = particle_set(iparticle)%r
          DO i = 1, 3
             particle_set(iparticle)%r(i) = rvec(i) + dx
             CALL get_ddapc(qs_env,.FALSE.,density_fit_section, qout1=qtot1,&
                  ext_rho_tot_g=rho_tot_g,Itype_of_density=type_of_density,error=error)
             particle_set(iparticle)%r(i) = rvec(i) - dx
             CALL get_ddapc(qs_env,.FALSE.,density_fit_section, qout1=qtot2,&
                  ext_rho_tot_g=rho_tot_g,Itype_of_density=type_of_density,error=error)
             ddqv = (qtot1 - qtot2)/(2.0_dp*dx)
             DO kk = 1, SIZE(qtot1)-1, SIZE(radii)
                IF (ANY(ddqv(kk:kk+2).GT.1.0E-8_dp)) THEN
                   WRITE(*,'(A,2F12.6,F12.2)')"Error :",SUM(dqv(kk:kk+2,iparticle,i)),SUM(ddqv(kk:kk+2)),&
                     ABS((SUM(ddqv(kk:kk+2))-SUM(dqv(kk:kk+2,iparticle,i)))/SUM(ddqv(kk:kk+2))*100.0_dp)
                END IF
             END DO
             particle_set(iparticle)%r = rvec
          END DO
       END DO
       !
       DEALLOCATE(qtot1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qtot2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ddqv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE debug_charge

END MODULE cp_ddapc_util
