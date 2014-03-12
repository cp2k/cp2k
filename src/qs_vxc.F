!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_vxc

  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE input_constants,                 ONLY: sic_ad,&
                                             sic_eo,&
                                             sic_mauri_spz,&
                                             sic_mauri_us,&
                                             sic_none,&
                                             xc_none,&
                                             xc_vdw_fun_nonloc
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grids,                        ONLY: pw_grid_compare
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE qs_dispersion_nonloc,            ONLY: calculate_dispersion_nonloc
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE xc,                              ONLY: xc_exc_calc,&
                                             xc_vxc_pw_create1
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: qs_vxc_create

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc'

CONTAINS

! *****************************************************************************
!> \brief calculates and allocates the xc potential, already reducing it to
!>      the dependence on rho and the one on tau
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param qs_env qs_environment to get all the needed things
!> \param rho_struct density for which v_xc is calculated
!> \param just_energy if true calculates just the energy, and does not
!>        allocate v_*_rspace
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      - 05.2002 modified to use the mp_allgather function each pe
!>        computes only part of the grid and this is broadcasted to all
!>        instead of summed.
!>        This scales significantly better (e.g. factor 3 on 12 cpus
!>        32 H2O) [Joost VdV]
!>      - moved to qs_ks_methods [fawzi]
!>      - sic alterations [Joost VandeVondele]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_vxc_create( qs_env, rho_struct, xc_section, vxc_rho, vxc_tau, exc,&
       just_energy, adiabatic_rescale_factor, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    REAL(KIND=dp)                            :: exc
    LOGICAL, INTENT(in), OPTIONAL            :: just_energy
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: adiabatic_rescale_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_vxc_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, myfun, stat, &
                                                vdw
    LOGICAL                                  :: do_adiabatic_rescaling, &
                                                failure, my_just_energy, &
                                                sic_scaling_b_zero, uf_grid, &
                                                vdW_nl
    REAL(KIND=dp)                            :: exc_m, factor, &
                                                my_adiabatic_rescale_factor, &
                                                my_scaling, nelec_s_inv, &
                                                nelec_spin(2)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_rho, my_vxc_tau, &
                                                rho_g, rho_m_gspace, &
                                                rho_m_rspace, rho_r, tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool, vdw_pw_pool, &
                                                xc_pw_pool
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_g2, tmp_pw
    TYPE(virial_type), POINTER               :: virial

    CALL timeset( routineN ,handle)

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
    NULLIFY(dft_control, pw_env, auxbas_pw_pool, xc_pw_pool, vdw_pw_pool, cell, my_vxc_rho, &
            tmp_pw, tmp_g, tmp_g2, my_vxc_tau, rho_g, rho_r, tau, rho_m_rspace, &
            rho_m_gspace)

    my_just_energy=.FALSE.
    IF (PRESENT(just_energy)) my_just_energy=just_energy
    my_adiabatic_rescale_factor = 1.0_dp
    do_adiabatic_rescaling = .FALSE.
    IF( PRESENT(adiabatic_rescale_factor)) THEN
       my_adiabatic_rescale_factor = adiabatic_rescale_factor
       do_adiabatic_rescaling = .TRUE.
    END IF

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
         pw_env=pw_env, cell=cell, virial=virial, error=error)

    CALL section_vals_val_get(xc_section,"XC_FUNCTIONAL%_SECTION_PARAMETERS_",&
         i_val=myfun,error=error)
    CALL section_vals_val_get(xc_section,"VDW_POTENTIAL%POTENTIAL_TYPE",&
         i_val=vdw,error=error)
    vdW_nl = (vdw==xc_vdw_fun_nonloc)

    IF (myfun/=xc_none .OR. vdW_nl) THEN

       ! test if the real space density is available
       CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_struct%rho_r_valid,cp_failure_level,routineP,error,failure)
       CALL cp_assert( dft_control%nspins == 1 .OR. dft_control%nspins == 2,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "nspins must be 1 or 2",error,failure)
       ! there are some options related to SIC here.
       ! Normal DFT computes E(rho_alpha,rho_beta) (or its variant E(2*rho_alpha) for non-LSD)
       ! SIC can             E(rho_alpha,rho_beta)-b*(E(rho_alpha,rho_beta)-E(rho_beta,rho_beta))
       ! or compute          E(rho_alpha,rho_beta)-b*E(rho_alpha-rho_beta,0)

       ! my_scaling is the scaling needed of the standard E(rho_alpha,rho_beta) term
       my_scaling=1.0_dp
       SELECT CASE (dft_control%sic_method_id)
       CASE ( sic_none )
          ! all fine
       CASE ( sic_mauri_spz, sic_ad )
          ! no idea yet what to do here in that case
          CPPrecondition(.NOT.rho_struct%tau_r_valid ,cp_failure_level,routineP,error,failure)
       CASE ( sic_mauri_us )
          my_scaling=1.0_dp-dft_control%sic_scaling_b
          ! no idea yet what to do here in that case
          CPPrecondition(.NOT.rho_struct%tau_r_valid ,cp_failure_level,routineP,error,failure)
       CASE ( sic_eo )
          ! NOTHING TO BE DONE
       CASE DEFAULT
          ! this case has not yet been treated here
          CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,"NYI",error,failure)
       END SELECT

       IF (dft_control%sic_scaling_b .EQ. 0.0_dp) THEN
          sic_scaling_b_zero = .TRUE.
       ELSE
          sic_scaling_b_zero = .FALSE.
       ENDIF

       IF ( .NOT. failure ) THEN
          CALL pw_env_get(pw_env,xc_pw_pool=xc_pw_pool,auxbas_pw_pool=auxbas_pw_pool,&
               error=error)
          uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid,xc_pw_pool%pw_grid)

          ALLOCATE(rho_r(dft_control%nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.not.uf_grid) THEN
             DO ispin=1,dft_control%nspins
                rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
             END DO

             IF (rho_struct%tau_r_valid) THEN
                ALLOCATE(tau(dft_control%nspins),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DO ispin=1,dft_control%nspins
                   tau(ispin)%pw => rho_struct%tau_r(ispin)%pw
                END DO
             END IF

             ! for gradient corrected functional the density in g space might
             ! be useful so if we have it, we pass it in
             IF ( rho_struct%rho_g_valid ) THEN
                ALLOCATE(rho_g(dft_control%nspins),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DO ispin=1,dft_control%nspins
                   rho_g(ispin)%pw => rho_struct%rho_g(ispin)%pw
                END DO
             END IF
          ELSE
             CPPrecondition(rho_struct%rho_g_valid,cp_failure_level,routineP,error,failure)
             ALLOCATE(rho_g(dft_control%nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,dft_control%nspins
                CALL pw_pool_create_pw(xc_pw_pool,rho_g(ispin)%pw,&
                     in_space=RECIPROCALSPACE, use_data=COMPLEXDATA1D,error=error)
                CALL pw_transfer(rho_struct%rho_g(ispin)%pw,rho_g(ispin)%pw, error=error)
             END DO
             DO ispin=1,dft_control%nspins
                CALL pw_pool_create_pw(xc_pw_pool,rho_r(ispin)%pw,&
                     in_space=REALSPACE, use_data=REALDATA3D,error=error)
                CALL pw_transfer(rho_g(ispin)%pw,rho_r(ispin)%pw, error=error)
             END DO
             IF (rho_struct%tau_r_valid) THEN
               ! tau with finer grids is not implemented (at least not correctly), which this asserts
               CALL cp_unimplemented_error(fromWhere=routineP, &
                    message="tau with finer grids", &
                    error=error, error_level=cp_failure_level)
!                ALLOCATE(tau(dft_control%nspins),stat=stat)
!                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!                DO ispin=1,dft_control%nspins
!                   CALL pw_pool_create_pw(xc_pw_pool,tau(ispin)%pw,&
!                        in_space=REALSPACE, use_data=REALDATA3D,error=error)
!
!                   CALL pw_pool_create_pw(xc_pw_pool,tmp_g,&
!                        in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
!                   CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g2,&
!                        in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
!                   CALL pw_transfer(tau(ispin)%pw,tmp_g, error=error)
!                   CALL pw_transfer(tmp_g,tmp_g2, error=error)
!                   CALL pw_transfer(tmp_g2,tmp_pw, error=error)
!                   CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g2,error=error)
!                   CALL pw_pool_give_back_pw(xc_pw_pool,tmp_g,error=error)
!                END DO
             END IF
          END IF

          ! add the nlcc densities
          IF (ASSOCIATED(qs_env%rho_nlcc)) THEN
             factor=1.0_dp
             DO ispin=1,dft_control%nspins
                CALL pw_axpy ( qs_env%rho_nlcc%pw, rho_r(ispin)%pw, factor, error)
                CALL pw_axpy ( qs_env%rho_nlcc_g%pw, rho_g(ispin)%pw, factor, error)
             ENDDO
          ENDIF

          !
          ! here the rho_r, rho_g, tau is what it should be
          ! we get back the right my_vxc_rho and my_vxc_tau as required
          !
          IF (my_just_energy) THEN
             exc=xc_exc_calc(rho_r=rho_r,tau=tau,&
                  rho_g=rho_g, xc_section=xc_section,&
                  cell=cell, pw_pool=xc_pw_pool,&
                  error=error)

          ELSE
             CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                  rho_g=rho_g,tau=tau,exc=exc,&
                  xc_section=xc_section,&
                  cell=cell, pw_pool=xc_pw_pool,&
                  virial=virial,&
                  error=error)
          END IF

          ! remove the nlcc densities (keep stuff in original state)
          IF (ASSOCIATED(qs_env%rho_nlcc)) THEN
             factor=-1.0_dp
             DO ispin=1,dft_control%nspins
                CALL pw_axpy ( qs_env%rho_nlcc%pw, rho_r(ispin)%pw, factor, error)
                CALL pw_axpy ( qs_env%rho_nlcc_g%pw, rho_g(ispin)%pw, factor, error)
             ENDDO
          ENDIF

          ! calclulate non-local vdW functional
          ! only if this XC_SECTION has it
          ! if yes, we use the dispersion_env from qs_env
          ! this is dangerous, as it assumes a special connection xc_section -> qs_env
          IF (vdW_nl) THEN
             ! no spin polarization allowed
             CPPrecondition(dft_control%nspins==1,cp_failure_level,routineP,error,failure)
             CALL pw_env_get(pw_env,vdw_pw_pool=vdw_pw_pool,error=error)
             IF(my_just_energy) THEN
                CALL calculate_dispersion_nonloc(qs_env,my_vxc_rho,rho_r,rho_g,my_just_energy,&
                  cell,vdw_pw_pool,xc_pw_pool,error=error)
             ELSE
                CALL calculate_dispersion_nonloc(qs_env,my_vxc_rho,rho_r,rho_g,my_just_energy,&
                  cell,vdw_pw_pool,xc_pw_pool,virial=virial,error=error)
             END IF
          END IF

          !! Apply rescaling to the potential if requested
          IF(.NOT. my_just_energy) THEN
             IF(do_adiabatic_rescaling) THEN
                IF( ASSOCIATED(my_vxc_rho)) THEN
                   DO ispin=1,SIZE(my_vxc_rho)
                      my_vxc_rho(ispin)%pw%cr3d=my_vxc_rho(ispin)%pw%cr3d*my_adiabatic_rescale_factor
                   END DO
                END IF
             END IF
          END IF

          IF (my_scaling .NE. 1.0_dp) THEN
             exc=exc * my_scaling
             IF (ASSOCIATED(my_vxc_rho)) THEN
                DO ispin=1,SIZE(my_vxc_rho)
                   my_vxc_rho(ispin)%pw%cr3d=my_vxc_rho(ispin)%pw%cr3d*my_scaling
                ENDDO
             ENDIF
             IF (ASSOCIATED(my_vxc_tau)) THEN
                DO ispin=1,SIZE(my_vxc_tau)
                   my_vxc_tau(ispin)%pw%cr3d=my_vxc_tau(ispin)%pw%cr3d*my_scaling
                ENDDO
             ENDIF
          ENDIF

          ! we have pw data for the xc, qs_ks requests coeff structure, here we transfer
          ! pw -> coeff
          IF (ASSOCIATED(my_vxc_rho)) THEN
             ALLOCATE(vxc_rho(dft_control%nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,dft_control%nspins
                vxc_rho(ispin)%pw => my_vxc_rho(ispin)%pw
             END DO
             DEALLOCATE(my_vxc_rho,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(my_vxc_tau)) THEN
             ALLOCATE(vxc_tau(dft_control%nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,dft_control%nspins
                vxc_tau(ispin)%pw => my_vxc_tau(ispin)%pw
             END DO
             DEALLOCATE(my_vxc_tau,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF

          ! compute again the xc but now for Exc(m,o) and the opposite sign
          IF (dft_control%sic_method_id .EQ. sic_mauri_spz .AND. .NOT. sic_scaling_b_zero) THEN
             ALLOCATE(rho_m_rspace(2),rho_m_gspace(2))
             CALL pw_pool_create_pw(xc_pw_pool, rho_m_gspace(1)%pw,&
                  use_data = COMPLEXDATA1D,&
                  in_space = RECIPROCALSPACE, error=error)
             CALL pw_pool_create_pw(xc_pw_pool, rho_m_rspace(1)%pw,&
                  use_data = REALDATA3D,&
                  in_space = REALSPACE, error=error)
             CALL pw_copy(rho_struct%rho_r(1)%pw,rho_m_rspace(1)%pw, error=error)
             CALL pw_axpy(rho_struct%rho_r(2)%pw,rho_m_rspace(1)%pw,alpha=-1._dp, error=error)
             CALL pw_copy(rho_struct%rho_g(1)%pw,rho_m_gspace(1)%pw, error=error)
             CALL pw_axpy(rho_struct%rho_g(2)%pw,rho_m_gspace(1)%pw,alpha=-1._dp, error=error)
             ! bit sad, these will be just zero...
             CALL pw_pool_create_pw(xc_pw_pool, rho_m_gspace(2)%pw,&
                  use_data = COMPLEXDATA1D,&
                  in_space = RECIPROCALSPACE, error=error)
             CALL pw_pool_create_pw(xc_pw_pool, rho_m_rspace(2)%pw,&
                  use_data = REALDATA3D,&
                  in_space = REALSPACE, error=error)
             CALL pw_zero(rho_m_rspace(2)%pw, error=error)
             CALL pw_zero(rho_m_gspace(2)%pw, error=error)

             rho_g(1)%pw => rho_m_gspace(1)%pw
             rho_g(2)%pw => rho_m_gspace(2)%pw
             rho_r(1)%pw => rho_m_rspace(1)%pw
             rho_r(2)%pw => rho_m_rspace(2)%pw

             IF (my_just_energy) THEN
                exc_m=xc_exc_calc(rho_r=rho_r,tau=tau,&
                     rho_g=rho_g, xc_section=xc_section,&
                     cell=cell, pw_pool=xc_pw_pool,&
                     error=error)
             ELSE
                CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                     rho_g=rho_g,tau=tau,exc=exc_m,&
                     xc_section=xc_section,&
                     cell=cell, pw_pool=xc_pw_pool,&
                     error=error)
             END IF

             exc = exc - dft_control%sic_scaling_b * exc_m

             ! and take care of the potential only vxc_rho is taken into account
             IF (.NOT. my_just_energy) THEN
                vxc_rho(1)%pw%cr3d=vxc_rho(1)%pw%cr3d-dft_control%sic_scaling_b *&
                     my_vxc_rho(1)%pw%cr3d
                vxc_rho(2)%pw%cr3d=vxc_rho(2)%pw%cr3d+dft_control%sic_scaling_b *&
                     my_vxc_rho(1)%pw%cr3d ! 1=m
                CALL pw_release(my_vxc_rho(1)%pw,error=error)
                CALL pw_release(my_vxc_rho(2)%pw,error=error)
                DEALLOCATE(my_vxc_rho,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF

             DO ispin=1,2
                CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_rspace(ispin)%pw,&
                     error=error)
                CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_gspace(ispin)%pw,&
                     error=error)
             ENDDO
             DEALLOCATE(rho_m_rspace)
             DEALLOCATE(rho_m_gspace)

          ENDIF

          ! now we have - sum_s N_s * Exc(rho_s/N_s,0)
          IF ( dft_control%sic_method_id .EQ. sic_ad  .AND. .NOT. sic_scaling_b_zero ) THEN

             ! find out how many elecs we have
             CALL get_qs_env(qs_env,mos=mo_array,error=error)
             CALL get_mo_set(mo_set=mo_array(1)%mo_set,n_el_f=nelec_spin(1))
             CALL get_mo_set(mo_set=mo_array(2)%mo_set,n_el_f=nelec_spin(2))

             ALLOCATE(rho_m_rspace(2),rho_m_gspace(2))
             DO ispin=1,2
                CALL pw_pool_create_pw(xc_pw_pool, rho_m_gspace(ispin)%pw,&
                     use_data = COMPLEXDATA1D,&
                     in_space = RECIPROCALSPACE, error=error)
                CALL pw_pool_create_pw(xc_pw_pool, rho_m_rspace(ispin)%pw,&
                     use_data = REALDATA3D,&
                     in_space = REALSPACE, error=error)
             ENDDO

             rho_g(1)%pw => rho_m_gspace(1)%pw
             rho_g(2)%pw => rho_m_gspace(2)%pw
             rho_r(1)%pw => rho_m_rspace(1)%pw
             rho_r(2)%pw => rho_m_rspace(2)%pw

             DO ispin=1,2
                IF (nelec_spin(ispin) .GT. 0.0_dp ) THEN
                   nelec_s_inv=1.0_dp/nelec_spin(ispin)
                ELSE
                   ! does it matter if there are no electrons with this spin (H) ?
                   nelec_s_inv=0.0_dp
                ENDIF
                CALL pw_copy(rho_struct%rho_r(ispin)%pw,rho_m_rspace(1)%pw, error=error)
                CALL pw_copy(rho_struct%rho_g(ispin)%pw,rho_m_gspace(1)%pw, error=error)
                CALL pw_scale(rho_m_rspace(1)%pw,nelec_s_inv, error=error)
                CALL pw_scale(rho_m_gspace(1)%pw,nelec_s_inv, error=error)
                CALL pw_zero(rho_m_rspace(2)%pw, error=error)
                CALL pw_zero(rho_m_gspace(2)%pw, error=error)

                IF (my_just_energy) THEN
                   exc_m=xc_exc_calc(rho_r=rho_r,tau=tau,&
                        rho_g=rho_g, xc_section=xc_section,&
                        cell=cell, pw_pool=xc_pw_pool,&
                        error=error)
                ELSE
                   CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                        rho_g=rho_g,tau=tau,exc=exc_m,&
                        xc_section=xc_section,&
                        cell=cell, pw_pool=xc_pw_pool,&
                        error=error)
                END IF

                exc = exc - dft_control%sic_scaling_b * nelec_spin(ispin) * exc_m

                ! and take care of the potential only vxc_rho is taken into account
                IF (.NOT. my_just_energy) THEN
                   vxc_rho(ispin)%pw%cr3d=vxc_rho(ispin)%pw%cr3d-dft_control%sic_scaling_b *&
                        my_vxc_rho(1)%pw%cr3d
                   CALL pw_release(my_vxc_rho(1)%pw,error=error)
                   CALL pw_release(my_vxc_rho(2)%pw,error=error)
                   DEALLOCATE(my_vxc_rho,stat=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ENDIF
             ENDDO

             DO ispin=1,2
                CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_rspace(ispin)%pw,&
                     error=error)
                CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_gspace(ispin)%pw,&
                     error=error)
             ENDDO
             DEALLOCATE(rho_m_rspace)
             DEALLOCATE(rho_m_gspace)

          ENDIF

          ! compute again the xc but now for Exc(n_down,n_down)
          IF (dft_control%sic_method_id .EQ. sic_mauri_us .AND. .NOT. sic_scaling_b_zero ) THEN
             rho_r(1)%pw => rho_struct%rho_r(2)%pw
             rho_r(2)%pw => rho_struct%rho_r(2)%pw
             IF ( rho_struct%rho_g_valid ) THEN
                rho_g(1)%pw => rho_struct%rho_g(2)%pw
                rho_g(2)%pw => rho_struct%rho_g(2)%pw
             ENDIF

             IF (my_just_energy) THEN
                exc_m=xc_exc_calc(rho_r=rho_r,tau=tau,&
                     rho_g=rho_g, xc_section=xc_section,&
                     cell=cell, pw_pool=xc_pw_pool,&
                     error=error)
             ELSE
                CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                     rho_g=rho_g,tau=tau,exc=exc_m,&
                     xc_section=xc_section,&
                     cell=cell, pw_pool=xc_pw_pool,&
                     error=error)
             END IF

             exc = exc + dft_control%sic_scaling_b * exc_m

             ! and take care of the potential
             IF (.NOT. my_just_energy) THEN
                ! both go to minority spin
                vxc_rho(2)%pw%cr3d = vxc_rho(2)%pw%cr3d + &
                     2.0_dp * dft_control%sic_scaling_b * my_vxc_rho(1)%pw%cr3d
                CALL pw_release(my_vxc_rho(1)%pw,error=error)
                CALL pw_release(my_vxc_rho(2)%pw,error=error)
                DEALLOCATE(my_vxc_rho)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF

          ENDIF

          !
          ! cleanups
          !
          IF (uf_grid) THEN
             DO ispin=1,SIZE(rho_r)
                CALL pw_pool_give_back_pw(xc_pw_pool,rho_r(ispin)%pw,error=error)
             END DO
             IF (ASSOCIATED(vxc_rho)) THEN
                DO ispin=1,SIZE(vxc_rho)
                   CALL pw_pool_create_pw(auxbas_pw_pool,tmp_pw,&
                        in_space=REALSPACE,use_data=REALDATA3D,error=error)

                   CALL pw_pool_create_pw(xc_pw_pool,tmp_g,&
                        in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                   CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g2,&
                        in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                   CALL pw_transfer(vxc_rho(ispin)%pw,tmp_g, error=error)
                   CALL pw_transfer(tmp_g,tmp_g2, error=error)
                   CALL pw_transfer(tmp_g2,tmp_pw, error=error)
                   CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g2,error=error)
                   CALL pw_pool_give_back_pw(xc_pw_pool,tmp_g,error=error)
                   !FM              CALL pw_zero(tmp_pw,error=error)
                   !FM              CALL pw_restrict_s3(vxc_rho(ispin)%pw,tmp_pw,&
                   !FM                   auxbas_pw_pool,param_section=interp_section,error=error)
                   CALL pw_pool_give_back_pw(xc_pw_pool,vxc_rho(ispin)%pw,error=error)
                   vxc_rho(ispin)%pw => tmp_pw
                   NULLIFY(tmp_pw)
                END DO
             END IF
             IF (ASSOCIATED(vxc_tau)) THEN
                DO ispin=1,SIZE(vxc_tau)
                   CALL pw_pool_create_pw(auxbas_pw_pool,tmp_pw,&
                        in_space=REALSPACE,use_data=REALDATA3D,error=error)

                   CALL pw_pool_create_pw(xc_pw_pool,tmp_g,&
                        in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                   CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g2,&
                        in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                   CALL pw_transfer(vxc_tau(ispin)%pw,tmp_g, error=error)
                   CALL pw_transfer(tmp_g,tmp_g2, error=error)
                   CALL pw_transfer(tmp_g2,tmp_pw, error=error)
                   CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g2,error=error)
                   CALL pw_pool_give_back_pw(xc_pw_pool,tmp_g,error=error)
                   !FM              CALL pw_zero(tmp_pw,error=error)
                   !FM              CALL pw_restrict_s3(vxc_rho(ispin)%pw,tmp_pw,&
                   !FM                   auxbas_pw_pool,param_section=interp_section,error=error)
                   CALL pw_pool_give_back_pw(xc_pw_pool,vxc_tau(ispin)%pw,error=error)
                   vxc_tau(ispin)%pw => tmp_pw
                   NULLIFY(tmp_pw)
                END DO
             END IF

          END IF
          DEALLOCATE(rho_r,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          IF (ASSOCIATED(rho_g)) THEN
             IF (uf_grid) THEN
                DO ispin=1,SIZE(rho_g)
                   CALL pw_pool_give_back_pw(xc_pw_pool,rho_g(ispin)%pw,error=error)
                END DO
             END IF
             DEALLOCATE(rho_g,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(tau)) THEN
             IF (uf_grid) THEN
                DO ispin=1,SIZE(tau)
                   CALL pw_pool_give_back_pw(xc_pw_pool,tau(ispin)%pw,error=error)
                END DO
             END IF
             DEALLOCATE(tau,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF

       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_vxc_create

END MODULE qs_vxc
