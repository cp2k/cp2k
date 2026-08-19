!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_twopt_thermodynamics
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   REAL(dp), PARAMETER, PUBLIC :: gas_constant = 8.31446261815324_dp
   REAL(dp), PARAMETER :: avogadro = 6.02214076E23_dp
   REAL(dp), PARAMETER :: boltzmann = 1.380649E-23_dp
   REAL(dp), PARAMETER :: light_speed_cm_s = 2.99792458E10_dp
   REAL(dp), PARAMETER :: planck = 6.62607015E-34_dp

   TYPE, PUBLIC :: twopt_channel_type
      REAL(dp) :: dof = 0.0_dp
      REAL(dp) :: gas_dof = 0.0_dp
      REAL(dp) :: s0_cm = 0.0_dp
      REAL(dp) :: k_parameter = 0.0_dp
      REAL(dp) :: fluidicity = 0.0_dp
      REAL(dp) :: packing_fraction = 0.0_dp
      REAL(dp) :: diffusion_cm2_s = 0.0_dp
      REAL(dp) :: zpe_kj_mol = 0.0_dp
      REAL(dp) :: energy_quantum_kj_mol = 0.0_dp
      REAL(dp) :: energy_classical_kj_mol = 0.0_dp
      REAL(dp) :: entropy_quantum_j_mol_k = 0.0_dp
      REAL(dp) :: entropy_classical_j_mol_k = 0.0_dp
      REAL(dp) :: free_energy_quantum_kj_mol = 0.0_dp
      REAL(dp) :: free_energy_classical_kj_mol = 0.0_dp
      REAL(dp) :: cv_quantum_j_mol_k = 0.0_dp
      REAL(dp) :: cv_classical_j_mol_k = 0.0_dp
   END TYPE twopt_channel_type

   TYPE, PUBLIC :: twopt_weights_type
      REAL(dp) :: entropy_translation = 0.0_dp
      REAL(dp) :: entropy_rotation = 0.0_dp
   END TYPE twopt_weights_type

   PUBLIC :: partition_twopt, build_twopt_weights, integrate_twopt_channel
   PUBLIC :: fluidicity_from_k, rotational_temperatures

CONTAINS

   REAL(dp) FUNCTION fluidicity_from_k(k_parameter)
      REAL(dp), INTENT(IN)                               :: k_parameter

      INTEGER                                            :: iteration
      REAL(dp)                                           :: derivative, f, f_new, polynomial, safe_k

      IF (k_parameter <= 0.0_dp) THEN
         fluidicity_from_k = 0.0_dp
         RETURN
      END IF
      safe_k = MAX(k_parameter, 1.0E-30_dp)
      f = MIN(MAX(0.7293_dp*safe_k**0.5727_dp, 1.0E-12_dp), 1.0_dp - 1.0E-12_dp)
      DO iteration = 1, 999
         polynomial = 2.0_dp*safe_k**(-4.5_dp)*f**7.5_dp - &
                      6.0_dp*safe_k**(-3.0_dp)*f**5.0_dp - &
                      safe_k**(-1.5_dp)*f**3.5_dp + &
                      6.0_dp*safe_k**(-1.5_dp)*f**2.5_dp + 2.0_dp*f - 2.0_dp
         derivative = 15.0_dp*safe_k**(-4.5_dp)*f**6.5_dp - &
                      30.0_dp*safe_k**(-3.0_dp)*f**4.0_dp - &
                      3.5_dp*safe_k**(-1.5_dp)*f**2.5_dp + &
                      15.0_dp*safe_k**(-1.5_dp)*f**1.5_dp + 2.0_dp
         IF (.NOT. ieee_is_finite(polynomial) .OR. .NOT. ieee_is_finite(derivative) .OR. &
             ABS(derivative) <= TINY(1.0_dp)) EXIT
         f_new = f - polynomial/derivative
         IF (.NOT. ieee_is_finite(f_new) .OR. f_new <= 0.0_dp .OR. f_new >= 1.0_dp) THEN
            IF (.NOT. ieee_is_finite(f_new) .OR. f_new <= 0.0_dp) THEN
               f_new = 0.5_dp*(f + 1.0E-12_dp)
            ELSE
               f_new = 0.5_dp*(f + 1.0_dp - 1.0E-12_dp)
            END IF
         END IF
         IF (ABS(f_new - f) <= 1.0E-12_dp) THEN
            f = f_new
            EXIT
         END IF
         f = f_new
      END DO
      fluidicity_from_k = MIN(MAX(f, 0.0_dp), 1.0_dp)
   END FUNCTION fluidicity_from_k

   SUBROUTINE partition_twopt(dos, frequency_step, temperature, molecular_mass, molecules, volume, &
                              gas, solid, result)
      REAL(dp), INTENT(IN)                               :: dos(0:), frequency_step, temperature, &
                                                            molecular_mass, volume
      INTEGER, INTENT(IN)                                :: molecules
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: gas(:), solid(:)
      TYPE(twopt_channel_type), INTENT(OUT)              :: result

      INTEGER                                            :: frequency
      REAL(dp)                                           :: gas_value, reference_packing, total_dof

      ALLOCATE (gas(0:UBOUND(dos, 1)), solid(0:UBOUND(dos, 1)))
      gas = 0.0_dp
      solid = dos
      result%s0_cm = MAX(dos(0), 0.0_dp)
      total_dof = trapezoid(dos, frequency_step)
      result%dof = total_dof
      IF (result%s0_cm <= 0.0_dp .OR. total_dof <= 0.0_dp .OR. molecules <= 0) RETURN

      result%k_parameter = (result%s0_cm/light_speed_cm_s*1.0E-2_dp)/REAL(molecules, dp)* &
                           SQRT(ACOS(-1.0_dp)*avogadro*boltzmann*temperature/ &
                                (molecular_mass*1.0E-3_dp))*2.0_dp/9.0_dp* &
                           (REAL(molecules, dp)/volume)**(1.0_dp/3.0_dp)*1.0E10_dp* &
                           (6.0_dp/ACOS(-1.0_dp))**(2.0_dp/3.0_dp)
      result%fluidicity = fluidicity_from_k(result%k_parameter)
      IF (result%fluidicity <= 0.0_dp) RETURN

      DO frequency = 0, UBOUND(dos, 1)
         gas_value = result%s0_cm/(1.0_dp + &
                     (ACOS(-1.0_dp)*result%s0_cm*frequency_step*REAL(frequency, dp)/ &
                      (6.0_dp*REAL(molecules, dp)*result%fluidicity))**2)
         gas(frequency) = MIN(MAX(gas_value, 0.0_dp), dos(frequency))
      END DO
      solid = MAX(dos - gas, 0.0_dp)
      result%gas_dof = trapezoid(gas, frequency_step)
      reference_packing = (result%fluidicity/result%k_parameter)**1.5_dp
      result%packing_fraction = reference_packing*result%gas_dof/total_dof
      IF (result%packing_fraction > 0.74_dp) THEN
         result%fluidicity = 0.0_dp
         result%packing_fraction = 0.0_dp
         result%gas_dof = 0.0_dp
         gas = 0.0_dp
         solid = dos
      END IF
   END SUBROUTINE partition_twopt

   SUBROUTINE build_twopt_weights(packing_fraction, molecular_mass, gas_molecules, temperature, volume, &
                                  rotational_temperature, rotational_symmetry, entropy_convention, weights)
      REAL(dp), INTENT(IN)                               :: packing_fraction, molecular_mass, gas_molecules, &
                                                            temperature, volume, rotational_temperature(3), &
                                                            rotational_symmetry
      CHARACTER(LEN=*), INTENT(IN)                       :: entropy_convention
      TYPE(twopt_weights_type), INTENT(OUT)              :: weights

      REAL(dp)                                           :: excess_entropy, inverse_lambda_cubed, mass_kg, &
                                                            packing, partition_rotation, z_factor

      weights%entropy_translation = 0.0_dp
      weights%entropy_rotation = 0.0_dp
      IF (gas_molecules <= 0.0_dp) RETURN
      packing = MIN(MAX(packing_fraction, 0.0_dp), 0.74_dp)
      mass_kg = molecular_mass*1.0E-3_dp/avogadro
      inverse_lambda_cubed = (2.0_dp*ACOS(-1.0_dp)*mass_kg*boltzmann*temperature/planck**2)**1.5_dp
      excess_entropy = packing*(3.0_dp*packing - 4.0_dp)/(1.0_dp - packing)**2
      IF (TRIM(entropy_convention) == "lin2003") THEN
         z_factor = (1.0_dp + packing + packing**2 - packing**3)/(1.0_dp - packing)**3
         excess_entropy = excess_entropy + LOG(z_factor)
      END IF
      weights%entropy_translation = (2.5_dp + &
         LOG(inverse_lambda_cubed*volume*1.0E-30_dp/gas_molecules) + excess_entropy)/3.0_dp

      IF (rotational_temperature(1) <= 0.0_dp) RETURN
      IF (rotational_temperature(3) < 0.0_dp) THEN
         partition_rotation = temperature/SQRT(rotational_temperature(1)*rotational_temperature(2))/ &
                              rotational_symmetry
         weights%entropy_rotation = (1.0_dp + LOG(partition_rotation))/2.0_dp
      ELSE
         partition_rotation = SQRT(ACOS(-1.0_dp)*temperature**3/PRODUCT(rotational_temperature))/ &
                              rotational_symmetry
         weights%entropy_rotation = (1.5_dp + LOG(partition_rotation))/3.0_dp
      END IF
   END SUBROUTINE build_twopt_weights

   SUBROUTINE integrate_twopt_channel(frequency, frequency_step, temperature, gas, solid, &
                                      gas_entropy_weight, result)
      REAL(dp), INTENT(IN)                               :: frequency(0:), frequency_step, temperature, &
                                                            gas(0:), solid(0:), gas_entropy_weight
      TYPE(twopt_channel_type), INTENT(INOUT)             :: result

      INTEGER                                            :: index
      REAL(dp)                                           :: classical_entropy_weight, classical_free_weight, &
                                                            cv_quantum_weight, energy_quantum_weight, &
                                                            free_quantum_weight, quantum_entropy_weight, &
                                                            weight

      DO index = 0, UBOUND(frequency, 1)
         CALL harmonic_weights(frequency(index), temperature, energy_quantum_weight, &
                               quantum_entropy_weight, free_quantum_weight, cv_quantum_weight, &
                               classical_entropy_weight, classical_free_weight)
         weight = integration_weight(index, UBOUND(frequency, 1))*frequency_step
         result%zpe_kj_mol = result%zpe_kj_mol + &
            weight*solid(index)*dimensionless_frequency(frequency(index), temperature)*0.5_dp* &
            temperature*gas_constant*1.0E-3_dp
         result%energy_quantum_kj_mol = result%energy_quantum_kj_mol + &
            weight*(solid(index)*energy_quantum_weight + 0.5_dp*gas(index))* &
            temperature*gas_constant*1.0E-3_dp
         result%energy_classical_kj_mol = result%energy_classical_kj_mol + &
            weight*(solid(index) + 0.5_dp*gas(index))*temperature*gas_constant*1.0E-3_dp
         result%entropy_quantum_j_mol_k = result%entropy_quantum_j_mol_k + &
            weight*(solid(index)*quantum_entropy_weight + gas(index)*gas_entropy_weight)*gas_constant
         result%entropy_classical_j_mol_k = result%entropy_classical_j_mol_k + &
            weight*(solid(index)*classical_entropy_weight + gas(index)*gas_entropy_weight)*gas_constant
         result%free_energy_quantum_kj_mol = result%free_energy_quantum_kj_mol + &
            weight*(solid(index)*free_quantum_weight + gas(index)*(0.5_dp - gas_entropy_weight))* &
            temperature*gas_constant*1.0E-3_dp
         result%free_energy_classical_kj_mol = result%free_energy_classical_kj_mol + &
            weight*(solid(index)*classical_free_weight + gas(index)*(0.5_dp - gas_entropy_weight))* &
            temperature*gas_constant*1.0E-3_dp
         result%cv_quantum_j_mol_k = result%cv_quantum_j_mol_k + &
            weight*(solid(index)*cv_quantum_weight + 0.5_dp*gas(index))*gas_constant
         result%cv_classical_j_mol_k = result%cv_classical_j_mol_k + &
            weight*(solid(index) + 0.5_dp*gas(index))*gas_constant
      END DO
   END SUBROUTINE integrate_twopt_channel

   SUBROUTINE harmonic_weights(frequency_cm, temperature, energy, entropy, free_energy, heat_capacity, &
                               classical_entropy, classical_free_energy)
      REAL(dp), INTENT(IN)                               :: frequency_cm, temperature
      REAL(dp), INTENT(OUT)                              :: energy, entropy, free_energy, heat_capacity, &
                                                            classical_entropy, classical_free_energy

      REAL(dp)                                           :: exponential, u

      u = dimensionless_frequency(frequency_cm, temperature)
      IF (u <= 0.0_dp) THEN
         energy = 0.0_dp
         entropy = 0.0_dp
         free_energy = 0.0_dp
         heat_capacity = 0.0_dp
         classical_entropy = 0.0_dp
         classical_free_energy = 0.0_dp
      ELSE IF (u > 500.0_dp) THEN
         energy = 0.5_dp*u
         entropy = 0.0_dp
         free_energy = 0.5_dp*u
         heat_capacity = 0.0_dp
         classical_entropy = 1.0_dp - LOG(u)
         classical_free_energy = LOG(u)
      ELSE
         exponential = EXP(u)
         energy = 0.5_dp*u + u/(exponential - 1.0_dp)
         entropy = u/(exponential - 1.0_dp) - LOG(1.0_dp - EXP(-u))
         free_energy = 0.5_dp*u + LOG(1.0_dp - EXP(-u))
         heat_capacity = (u/(2.0_dp*SINH(0.5_dp*u)))**2
         classical_entropy = 1.0_dp - LOG(u)
         classical_free_energy = LOG(u)
      END IF
   END SUBROUTINE harmonic_weights

   PURE REAL(dp) FUNCTION dimensionless_frequency(frequency_cm, temperature)
      REAL(dp), INTENT(IN)                               :: frequency_cm, temperature

      dimensionless_frequency = planck*light_speed_cm_s*frequency_cm/(boltzmann*temperature)
   END FUNCTION dimensionless_frequency

   PURE REAL(dp) FUNCTION trapezoid(values, spacing)
      REAL(dp), INTENT(IN)                               :: values(0:), spacing

      INTEGER                                            :: last

      last = UBOUND(values, 1)
      IF (last == 0) THEN
         trapezoid = 0.0_dp
      ELSE
         trapezoid = spacing*(0.5_dp*(values(0) + values(last)) + SUM(values(1:last - 1)))
      END IF
   END FUNCTION trapezoid

   PURE REAL(dp) FUNCTION integration_weight(index, last)
      INTEGER, INTENT(IN)                                :: index, last

      IF (index == 0 .OR. index == last) THEN
         integration_weight = 0.5_dp
      ELSE
         integration_weight = 1.0_dp
      END IF
   END FUNCTION integration_weight

   SUBROUTINE rotational_temperatures(inertia_amu_angstrom2, linear, temperatures)
      REAL(dp), INTENT(IN)                               :: inertia_amu_angstrom2(3)
      LOGICAL, INTENT(IN)                                :: linear
      REAL(dp), INTENT(OUT)                              :: temperatures(3)

      INTEGER                                            :: axis, source
      REAL(dp), PARAMETER                                :: amu_angstrom2_to_kg_m2 = 1.66053906660E-47_dp

      temperatures = 0.0_dp
      IF (linear) THEN
         DO axis = 1, 2
            source = axis + 1
            IF (inertia_amu_angstrom2(source) > 0.0_dp) &
               temperatures(axis) = planck**2/(8.0_dp*ACOS(-1.0_dp)**2* &
                                    inertia_amu_angstrom2(source)*amu_angstrom2_to_kg_m2*boltzmann)
         END DO
         temperatures(3) = -1.0_dp
      ELSE
         DO axis = 1, 3
            IF (inertia_amu_angstrom2(axis) > 0.0_dp) &
               temperatures(axis) = planck**2/(8.0_dp*ACOS(-1.0_dp)**2* &
                                    inertia_amu_angstrom2(axis)*amu_angstrom2_to_kg_m2*boltzmann)
         END DO
      END IF
   END SUBROUTINE rotational_temperatures

END MODULE trajana_twopt_thermodynamics
