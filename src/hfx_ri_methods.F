!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the master routines for the Hartree-Fock RI Method
!> \author JGH July 2009
! *****************************************************************************
MODULE hfx_ri_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: geminal_basis_set_type,&
                                             get_geminal_basis_set,&
                                             write_geminal_basis_set
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_geminals,                     ONLY: geminal_charge,&
                                             geminal_coulomb,&
                                             geminal_coulomb_diaginv,&
                                             geminal_gto_coulomb
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE ri_environment_types,            ONLY: &
       get_ri_env, ri_env_create, ri_environment_type, ri_vector_add, &
       ri_vector_copy, ri_vector_create, ri_vector_dot, ri_vector_norm, &
       ri_vector_release, ri_vector_scale, ri_vector_type, set_ri_env
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_create,&
                                             semi_empirical_si_finalize,&
                                             semi_empirical_si_initialize,&
                                             semi_empirical_si_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_methods'

  PUBLIC :: hfx_ri_energy_potential, hfx_ri_env_create, hfx_ri_energy

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an ri_env
!> \param qs_env ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE hfx_ri_env_create (qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, istat, &
                                                maxlgem(2), nat, ngem, nkind, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natoms, nbasis
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff, ri_rhs
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(semi_empirical_si_type), POINTER    :: ri_store_int1, ri_store_int2, &
                                                ri_store_int3

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)

    CALL ri_env_create(ri_env,error)

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,error=error)
    nkind = SIZE(qs_kind_set)
    ALLOCATE(nbasis(nkind),natoms(nkind),stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    nbasis=0
    natoms=0
    DO ikind=1,nkind
       NULLIFY(geminal_basis_set)
       CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
       CALL get_qs_kind(qs_kind_set(ikind), geminal_basis_set=geminal_basis_set)
       IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=ngem)
       natoms(ikind) = nat
       nbasis(ikind) = ngem
    END DO

    NULLIFY(ri_coeff,ri_rhs)
    CALL ri_vector_create(ri_coeff, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_rhs, nkind, natoms, nbasis, error)

    CALL set_ri_env(ri_env,coeff=ri_coeff,rhs=ri_rhs,error=error)
    CALL set_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,error=error)

    !integral storage
    NULLIFY(ri_store_int1,ri_store_int2,ri_store_int3)
    hfx_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%HF",error=error)
    CALL semi_empirical_si_create(ri_store_int1, hfx_section, compression=.TRUE., error=error)
    CALL set_ri_env(ri_env, gem_integral_storage=ri_store_int1, error=error)
    CALL semi_empirical_si_create(ri_store_int2, hfx_section, compression=.TRUE., error=error)
    CALL set_ri_env(ri_env, ggg_integral_storage=ri_store_int2, error=error)
    CALL semi_empirical_si_create(ri_store_int3, hfx_section, compression=.FALSE., error=error)
    CALL set_ri_env(ri_env, gemdiag_integral_storage=ri_store_int3, error=error)

    !fraction of exchange
    CALL section_vals_val_get(hfx_section, "FRACTION", r_val=ri_env%xfraction, error=error)

    !input
    hfx_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%HF%HFX_RI",error=error)

    CALL section_vals_val_get(hfx_section,"EPS_SCREENING",r_val=ri_env%screening,error=error)
    CALL section_vals_val_get(hfx_section,"EPS_OPTIMIZATION",r_val=ri_env%eps_opt,error=error)
    CALL section_vals_val_get(hfx_section,"MAX_ITER",i_val=ri_env%max_iter,error=error)

    output_unit=cp_print_key_unit_nr(logger,hfx_section,"PRINT%GEMINAL_BASIS",&
         extension=".RILog",error=error)
    IF ( output_unit > 0 ) THEN
       DO ikind=1,nkind
          NULLIFY(geminal_basis_set)
          CALL get_qs_kind(qs_kind_set(ikind), geminal_basis_set=geminal_basis_set)
          IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
          CALL write_geminal_basis_set(geminal_basis_set,output_unit,error)
       END DO
    END IF

    output_unit = cp_logger_get_default_io_unit(logger)
    IF ( output_unit > 0 ) THEN
       CALL get_qs_kind_set(qs_kind_set=qs_kind_set,maxlgem=maxlgem)
       ngem=0
       DO ikind=1,nkind
         ngem = ngem + nbasis(ikind)*natoms(ikind)
       END DO
       WRITE(output_unit,'(/,A,A,T79,I2)') "  HFX_RI|",&
             " Maximum l quantum number for R variable",maxlgem(1)
       WRITE(output_unit,'(A,A,T79,I2)') "  HFX_RI|",&
             " Maximum l quantum number for S variable",maxlgem(2)
       WRITE(output_unit,'(A,A,T71,I10)') "  HFX_RI|",&
             " Total number of Geminal basis functions",ngem
    END IF

    DEALLOCATE(nbasis,natoms,stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    CALL timestop(handle)

  END SUBROUTINE hfx_ri_env_create

! *****************************************************************************
!> \brief Calculate energy and potential for the HFX RI approximation
!>
!> \param qs_env ...
!> \param calculate_force ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE hfx_ri_energy_potential (qs_env,calculate_force,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_energy_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: calculate_energy, &
                                                calculate_fock, &
                                                calculate_virial, &
                                                s_mstruct_changed
    REAL(KIND=dp)                            :: t1, t2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    t1 = m_walltime()
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)

    NULLIFY(ks_env)
    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    s_mstruct_changed=s_mstruct_changed,&
                    error=error)

    ! calcluate expansion coefficients
    CALL optimize_ri_coeff (qs_env,s_mstruct_changed,.TRUE.,error)

    CALL get_qs_env(qs_env=qs_env,virial=virial,error=error)
    calculate_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    ! calcluate Fock matrix and energy
    calculate_energy = .TRUE.
    calculate_fock = .TRUE.
    CALL ri_energy_potential(qs_env,calculate_energy,calculate_fock,calculate_force,calculate_virial,error)

    t2 = m_walltime()
    IF ( output_unit > 0 ) THEN
       WRITE(output_unit,'(T15,A,T25,A,T69,F12.3)') "HFX_RI|","Total time for calculation ",t2-t1
    END IF

    CALL timestop(handle)

  END SUBROUTINE hfx_ri_energy_potential

! *****************************************************************************
!> \brief calculates the HFX energy in RI approximation
!> \param qs_env ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE hfx_ri_energy (qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    ! calcluate expansion coefficients
    CALL optimize_ri_coeff (qs_env,.TRUE.,.FALSE.,error)

    ! calcluate energy
    CALL ri_energy_potential(qs_env,.TRUE.,.FALSE.,.FALSE.,.FALSE.,error)

    CALL timestop(handle)

  END SUBROUTINE hfx_ri_energy

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param s_mstruct_changed ...
!> \param print_on ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE optimize_ri_coeff (qs_env,s_mstruct_changed,print_on,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: s_mstruct_changed, print_on
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_ri_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, istat, iter, &
                                                niter, nkind, ounit, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natoms, nbasis
    LOGICAL                                  :: converged
    REAL(KIND=dp)                            :: alpha, beta, charge, energy, &
                                                eps_ri, papk, rnorm, rrk, &
                                                rro, t1, t2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_apk, ri_coeff, ri_pk, &
                                                ri_res, ri_rhs, ri_zk
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(semi_empirical_si_type), POINTER    :: gem_integral_storage, &
                                                gemdiag_integral_storage, &
                                                ggg_integral_storage

    CALL timeset(routineN,handle)

    t1 = m_walltime()
    logger => cp_error_get_logger(error)

    hfx_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%HF%HFX_RI",error=error)
    output_unit=cp_print_key_unit_nr(logger,hfx_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".RILog",error=error)
    IF ( output_unit > 0 .AND. print_on ) THEN
       WRITE(output_unit,'(A)') "  Hartree-Fock-Exchange Resolution of Identity Approximation"
    END IF

    NULLIFY(ri_env)
    CALL get_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,error=error)

    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,rhs=ri_rhs,error=error)

    ! integral storage
    CALL get_ri_env(ri_env=ri_env,gem_integral_storage=gem_integral_storage,&
                    ggg_integral_storage=ggg_integral_storage,&
                    gemdiag_integral_storage=gemdiag_integral_storage,error=error)

    ! calculate rhs of linear equation
    CALL semi_empirical_si_initialize(ggg_integral_storage,s_mstruct_changed,error)
    CALL geminal_gto_coulomb(qs_env,.FALSE.,.FALSE.,.FALSE.,.FALSE.,error)
    CALL semi_empirical_si_finalize(ggg_integral_storage,s_mstruct_changed,error)

    !allocate scratch vectors
    nkind = SIZE(ri_coeff%vector)
    ALLOCATE(nbasis(nkind),natoms(nkind),stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    DO ikind=1,nkind
      natoms(ikind) = ri_coeff%vector(ikind)%natom
      nbasis(ikind) = ri_coeff%vector(ikind)%nbasis
    END DO

    NULLIFY(ri_res,ri_zk,ri_pk,ri_apk)
    CALL ri_vector_create(ri_res, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_zk, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_pk, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_apk, nkind, natoms, nbasis, error)

    DEALLOCATE(nbasis,natoms,stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    niter = ri_env%max_iter
    eps_ri = ri_env%eps_opt

    converged = .FALSE.
    ! calculate residual
    CALL semi_empirical_si_initialize(gem_integral_storage,s_mstruct_changed,error)
    CALL geminal_coulomb(ri_coeff,ri_res,qs_env,energy,.FALSE.,.FALSE.,.FALSE.,error)
    CALL semi_empirical_si_finalize(gem_integral_storage,s_mstruct_changed,error)

    CALL ri_vector_add(-1._dp, ri_rhs, ri_res, error)
    CALL ri_vector_scale(-1._dp, ri_res, error)
    ! apply preconditioner
    CALL semi_empirical_si_initialize(gemdiag_integral_storage,s_mstruct_changed,error)
    CALL geminal_coulomb_diaginv(ri_res,ri_zk,qs_env,error)
    CALL semi_empirical_si_finalize(gemdiag_integral_storage,s_mstruct_changed,error)
    !
    CALL ri_vector_copy(ri_zk, ri_pk, error)
    CALL ri_vector_dot(rrk, ri_res, ri_zk, error)
    DO iter=1,niter
       CALL semi_empirical_si_initialize(gem_integral_storage,.FALSE.,error)
       CALL geminal_coulomb(ri_pk,ri_apk,qs_env,energy,.FALSE.,.FALSE.,.FALSE.,error)
       CALL semi_empirical_si_finalize(gem_integral_storage,.FALSE.,error)
       CALL ri_vector_dot(papk, ri_pk, ri_apk, error)
       alpha = rrk/papk
       CALL ri_vector_add(alpha, ri_pk, ri_coeff, error)
       CALL ri_vector_add(-alpha, ri_apk, ri_res, error)
       CALL ri_vector_norm(rnorm, ri_res, error)

       IF ( output_unit > 0 .AND. print_on ) THEN
          WRITE(*,'(T15,A,T25,A,T35,i4,T50,A,T61,F20.10)') "HFX_RI|","Iteration :",iter,"Residual=",rnorm
       END IF

       IF ( rnorm < eps_ri ) THEN
          converged = .TRUE.
          EXIT
       END IF
       ! apply preconditioner
       CALL semi_empirical_si_initialize(gemdiag_integral_storage,.FALSE.,error)
       CALL geminal_coulomb_diaginv(ri_res,ri_zk,qs_env,error)
       CALL semi_empirical_si_finalize(gemdiag_integral_storage,.FALSE.,error)
       !
       rro = rrk
       CALL ri_vector_dot(rrk, ri_res, ri_zk, error)
       beta = rrk/rro
       CALL ri_vector_scale(beta, ri_pk, error)
       CALL ri_vector_add(1._dp, ri_zk, ri_pk, error)
    END DO

    CALL ri_vector_release(ri_res, error)
    CALL ri_vector_release(ri_zk, error)
    CALL ri_vector_release(ri_pk, error)
    CALL ri_vector_release(ri_apk, error)

    ounit=cp_print_key_unit_nr(logger,hfx_section,"PRINT%CHARGE",&
         extension=".RILog",error=error)
    IF ( ounit > 0 .AND. print_on ) THEN
       charge=0._dp
       CALL geminal_charge(ri_coeff,qs_env,charge,error)
       WRITE(output_unit,'(T15,A,T25,A,T57,F24.12)') "HFX_RI|","Total charge of RI expansion ",charge
    END IF

    t2 = m_walltime()
    IF ( output_unit > 0 .AND. print_on ) THEN
       WRITE(output_unit,'(T15,A,T25,A,T69,F12.3)') "HFX_RI|","Total time for optimization ",t2-t1
    END IF

    CALL timestop(handle)

  END SUBROUTINE optimize_ri_coeff

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_energy ...
!> \param calculate_fock ...
!> \param calculate_force ...
!> \param calculate_virial ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ri_energy_potential (qs_env,calculate_energy,calculate_fock,&
                                  calculate_force,calculate_virial,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_energy, &
                                                calculate_fock, &
                                                calculate_force, &
                                                calculate_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_energy_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, istat, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natoms, nbasis
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff, ri_res, ri_rhs
    TYPE(semi_empirical_si_type), POINTER    :: gem_integral_storage, &
                                                ggg_integral_storage

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,energy=energy,error=error)
    ! integral storage
    CALL get_ri_env(ri_env=ri_env,gem_integral_storage=gem_integral_storage,&
                    ggg_integral_storage=ggg_integral_storage,error=error)

    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,rhs=ri_rhs,error=error)
    nkind = SIZE(ri_coeff%vector)
    ALLOCATE(nbasis(nkind),natoms(nkind),stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    DO ikind=1,nkind
      natoms(ikind) = ri_coeff%vector(ikind)%natom
      nbasis(ikind) = ri_coeff%vector(ikind)%nbasis
    END DO

    NULLIFY(ri_res)
    CALL ri_vector_create(ri_res, nkind, natoms, nbasis, error)

    DEALLOCATE(nbasis,natoms,stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    CALL semi_empirical_si_initialize(gem_integral_storage,.FALSE.,error)
    CALL geminal_coulomb(ri_coeff,ri_res,qs_env,ri_env%ehfx1,calculate_energy,calculate_force,&
                         calculate_virial,error)
    CALL semi_empirical_si_finalize(gem_integral_storage,.FALSE.,error)

    CALL semi_empirical_si_initialize(ggg_integral_storage,.FALSE.,error)
    CALL geminal_gto_coulomb(qs_env,calculate_fock,calculate_energy,&
                             calculate_force,calculate_virial,error)
    CALL semi_empirical_si_finalize(ggg_integral_storage,.FALSE.,error)

    energy%ex = 0.5_dp*(ri_env%ehfx1 + ri_env%ehfx2)

    CALL ri_vector_release(ri_res, error)

    CALL timestop(handle)

  END SUBROUTINE ri_energy_potential

! *****************************************************************************

END MODULE hfx_ri_methods

