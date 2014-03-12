!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains some functions used in the context of adiabatic hybrid functionals
!> \par History
!>      01.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE xc_adiabatic_methods
  USE f77_blas
  USE input_constants,                 ONLY: do_hfx_potential_coulomb
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: oorootpi
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC     rescale_MCY3_pade

  CONTAINS

! *****************************************************************************
!> \brief - Calculates rescaling factors for XC potentials and energy expression
!> \param qs_en v
!> \param hf_energy Array of size 2 containing the two HF energies (Ex^{HF} and Ex^{HF,LR}
!> \param energy QS energy type
!> \param adiabatic_lambda , adiabatic_omega: Parameters for adiabatic connection
!> \param scale_ * : scaling coefficients for xc-potentials to be calculated
!> \param total_energy will contain the full xc energy
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Model for adiabatic connection:
!>
!>         W_lambda = a + b*lambda/(1+c*lambda)
!>         Exc = a + b*(c-log(1+c)/c^2)
!>         a = Ex^{HF}
!>         b = -c1*2*omega/sqrt(PI)*nelectron
!>         c = -1/lambda - b/(Ex^{HF}-W_lambda^{BLYP} + c2*W_lambda^{B88,LR}-c3*W_lambda^{HF,LR}
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE rescale_MCY3_pade(qs_env, hf_energy, energy, adiabatic_lambda, &
                               adiabatic_omega, scale_dEx1, scale_ddW0, scale_dDFA,&
                               scale_dEx2, total_energy_xc, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp), INTENT(INOUT)                  :: hf_energy(*)
    TYPE(qs_energy_type), POINTER            :: energy
    REAL(dp), INTENT(IN)                     :: adiabatic_lambda, &
                                                adiabatic_omega
    REAL(dp), INTENT(INOUT)                  :: scale_dEx1, scale_ddW0, &
                                                scale_dDFA, scale_dEx2, &
                                                total_energy_xc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: nelec_a, nelec_b, nelectron, &
                                                nspins
    LOGICAL                                  :: do_swap_hf
    REAL(dp) :: a, b, c, c1, da_dDFA, da_ddW0, da_dEx1, da_dEx2, db_dDFA, &
      db_ddW0, db_dEx1, db_dEx2, dc_dDFA, dc_ddW0, dc_dEx1, dc_dEx2, dExc_da, &
      dExc_db, dExc_dc, dfa_energy, swap_value
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    do_swap_hf = .FALSE.
    !! Assume the first HF section is the Coulomb one
    IF( qs_env%x_data(1,1)%potential_parameter%potential_type /= do_hfx_potential_coulomb ) do_swap_hf = .TRUE.

    IF( do_swap_hf ) THEN
      swap_value = hf_energy(1)
      hf_energy(1) = hf_energy(2)
      hf_energy(2) = swap_value
    END IF

    c1 = 0.23163_dp

    CALL get_qs_env(qs_env=qs_env, mos=mos, error=error)
    CALL get_mo_set(mo_set=mos(1)%mo_set,nelectron=nelec_a)
    nspins = SIZE(mos)
    IF(nspins == 2 ) THEN
      CALL get_mo_set(mo_set=mos(2)%mo_set,nelectron=nelec_b)
    ELSE
      nelec_b = 0
    END IF
    nelectron = nelec_a + nelec_b
    dfa_energy = energy%exc + energy%exc1
    a = hf_energy(1)
    b = -c1*2.0_dp*adiabatic_omega*oorootpi*nelectron!-0.23163_dp*2.0_dp*0.2_dp*oorootpi*nelectron
    c = -1.0_dp / adiabatic_lambda -b/(hf_energy(1)-dfa_energy-hf_energy(2))

    dExc_da = 1.0_dp
    dExc_db = 1.0_dp/c -(LOG(ABS(1.0_dp + c))/(c*c))
    dExc_dc = -b/(c*c*c*(1.0_dp+c))*(2.0_dp*c+c*c-2.0_dp*LOG(ABS(1.0_dp+c))-2.0_dp*LOG(ABS(1.0_dp+c))*c)

    da_dEx1 = 1.0_dp
    da_ddW0 = 0.0_dp
    da_dDFA = 0.0_dp
    da_dEx2 = 0.0_dp

    db_dEx1 = 0.0_dp
    db_ddW0 = 1.0_dp
    db_dDFA = 0.0_dp
    db_dEx2 = 0.0_dp

    dc_dEx1 =  b / (hf_energy(1)-dfa_energy-hf_energy(2))**2
    dc_ddW0 = -1.0_dp / (hf_energy(1)-dfa_energy-hf_energy(2))
    dc_dDFA = -dc_dEx1
    dc_dEx2 = -dc_dEx1

    scale_dEx1 = dExc_da*da_dEx1 + dExc_db*db_dEx1 + dExc_dc*dc_dEx1
    scale_ddW0 = dExc_da*da_ddW0 + dExc_db*db_ddW0 + dExc_dc*dc_ddW0
    scale_dDFA = dExc_da*da_dDFA + dExc_db*db_dDFA + dExc_dc*dc_dDFA
    scale_dEx2 = dExc_da*da_dEx2 + dExc_db*db_dEx2 + dExc_dc*dc_dEx2

    total_energy_xc = a + b/(c*c)*(c-LOG(ABS(1.0_dp+c)))
    IF( do_swap_hf ) THEN
      swap_value = scale_dEx1
      scale_dEx1 = scale_dEx2
      scale_dEx2 = swap_value
    END IF
  END SUBROUTINE rescale_MCY3_pade

END MODULE xc_adiabatic_methods
