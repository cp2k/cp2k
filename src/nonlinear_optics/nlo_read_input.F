!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \author Jan Wilhelm, with help from ChatGPT

MODULE nlo_read_input
   USE bibliography,                    ONLY: Cheng2015
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              real_t,&
                                              val_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE nlo_types,                       ONLY: nlo_frequency_type,&
                                              nlo_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo_read_input'

   PUBLIC :: create_nonlinear_optics_section, read_nlo_input
   PUBLIC :: nlo_frequency_case_count, nlo_frequency_case

CONTAINS

! **************************************************************************************************
!> \brief Define the nonlinear susceptibility model, conventions, and input parameters.
!> \param section Created NONLINEAR_OPTICS section of the input tree.
! **************************************************************************************************
   SUBROUTINE create_nonlinear_optics_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(LEN=6)                                   :: frequency_name
      INTEGER                                            :: i_frequency
      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="NONLINEAR_OPTICS", &
                          description= &
                          "Frequency-dependent nonlinear optical response of an insulating "// &
                          "periodic Kohn-Sham system in the independent-particle "// &
                          "approximation (no explicit excited-state treatment, in particular"// &
                          " no TDDFT, GW-BSE or excitonic effects). The induced polarization"// &
                          " is expanded in powers of the electric field. With $E(t)=\sum_ω "// &
                          "E(ω)e^{-iωt}"//"$, the second- through fifth-order "// &
                          "contributions are "// &
                          ACHAR(10)//ACHAR(10)//"$$P_i^{(2)}(Ω)=ε_0\sum_{ab}\sum_{ω_1+ω_2=Ω} "// &
                          "χ_{i;ab}^{(2)}(ω_1,ω_2)E_a(ω_1)E_b(ω_2),$$"//ACHAR(10)//ACHAR(10)// &
                          ACHAR(10)//ACHAR(10)// &
                          "$$P_i^{(3)}(Ω)=ε_0\sum_{abc} \sum_{ω_1+ω_2+ω_3=Ω} "// &
                          "χ_{i;abc}^{(3)}(ω_1,ω_2,ω_3) E_a(ω_1)E_b(ω_2)E_c(ω_3),$$"//ACHAR(10)// &
                          ACHAR(10)//ACHAR(10)//ACHAR(10)// &
                          "$$P_i^{(4)}(Ω)=ε_0\sum_{abcd} \sum_{ω_1+ω_2+ω_3+ω_4=Ω} "// &
                          "χ_{i;abcd}^{(4)}(ω_1,ω_2,ω_3,ω_4) "// &
                          "E_a(ω_1)E_b(ω_2)E_c(ω_3)E_d(ω_4),$$"//ACHAR(10)//ACHAR(10)//ACHAR(10)// &
                          ACHAR(10)// &
                          "$$P_i^{(5)}(Ω)=ε_0\sum_{abcde} \sum_{ω_1+ω_2+ω_3+ω_4+ω_5=Ω} "// &
                          "χ_{i;abcde}^{(5)}(ω_1,ω_2,ω_3,ω_4,ω_5) "// &
                          "E_a(ω_1)E_b(ω_2)E_c(ω_3)E_d(ω_4)E_e(ω_5).$$"//ACHAR(10)//ACHAR(10)// &
                          "Here i,a,... denote Cartesian directions and ε_0 is the vacuum "// &
                          "permittivity. Each tensor is averaged over all simultaneous "// &
                          "permutations of its incident frequency-direction pairs. The "// &
                          "frequency sums include ordered tuples, so no extra factorial "// &
                          "occurs in these formulas. Negative frequencies describe "// &
                          "complex-conjugate fields; harmonic generation uses equal positive"// &
                          " frequencies. The input specifies photon energies ℏω_j in eV, "// &
                          "with Ω≠0. The density matrix is expanded recursively in the "// &
                          "fields, using all AO bands and analytic derivatives with respect "// &
                          "to k. In the energy eigenbasis the usual "// &
                          "population and coherence relaxation terms are "//ACHAR(10)//ACHAR(10)// &
                          "$$\left.\dot{ρ}_{mn}\right|_{\mathrm{rel}}= "// &
                          "-δ_{mn}(ρ_{nn}-f_n)/T_1 -(1-δ_{mn})ρ_{mn}/T_2.$$"//ACHAR(10)// &
                          ACHAR(10)// &
                          "Here $f_n$ is the equilibrium occupation. This implementation "// &
                          "uses $T_1=T_2=T$, giving the basis-covariant "// &
                          "relaxation operator and common energy broadening "//ACHAR(10)// &
                          ACHAR(10)//"$$\left.\dot{ρ}\right|_{\mathrm{rel}}= "// &
                          "-(ρ-ρ_{\mathrm{eq}})/T, \qquad Γ=ℏ/T,\qquad T_2=ℏ/Γ.$$"//ACHAR(10)// &
                          ACHAR(10)// &
                          "At each response step, denominators are $ε_m-ε_n-ℏΩ-iΓ$. The "// &
                          "current includes the field-dependent contribution "// &
                          "required by the projected position operator. In "// &
                          "atomic units, with r_i=i∂_{k_i}+A_i, "//ACHAR(10)//ACHAR(10)// &
                          "$$F_{ia}=∂_{k_i}A_a-∂_{k_a}A_i-i[A_i,A_a], \qquad "// &
                          "v_i(E)=i[H+E_ar_a, r_i]=v_i(0)+E_aF_{ia}.$$"//ACHAR(10)//ACHAR(10)// &
                          "Let I^{(n)} denote the cell-integrated current response, "// &
                          "including the Brillouin-zone weights and spin "// &
                          "occupations, in SI units. The volume normalization and "// &
                          "conversion to susceptibility are "//ACHAR(10)//ACHAR(10)// &
                          "$$V=|\det(\mathbf a_1,\mathbf a_2,\mathbf a_3)|, \qquad "// &
                          "j^{(n)}=I^{(n)}/V,\qquad χ_{\mathrm{3D}}^{(n)}= "// &
                          "σ_{\mathrm{3D}}^{(n)}/(-iΩε_0),$$"//ACHAR(10)//ACHAR(10)// &
                          "where $σ^{(n)}"//"$ is the coefficient of the n field factors in "// &
                          "$j^{(n)}"//"$. For a two-dimensional slab with vacuum along "// &
                          "nonperiodic z, the following sheet quantities remove the "// &
                          "dependence on the chosen vacuum height: "//ACHAR(10)//ACHAR(10)// &
                          "$$A=|\mathbf a_1\times\mathbf a_2|,\qquad L_z=V/A, \qquad "// &
                          "P_{\mathrm{2D}}^{(n)}=L_zP_{\mathrm{3D}}^{(n)}, \qquad "// &
                          "χ_{\mathrm{2D}}^{(n)}=L_z χ_{\mathrm{3D}}^{(n)}.$$"//ACHAR(10)// &
                          ACHAR(10)// &
                          "The units are $(m/V)^{n-1}"//"$ for $χ_{3D}"//"$ and "// &
                          "$m^n/V^{n-1}"//"$ for "// &
                          "$χ_{2D}"//"$. Converge the basis set, the k mesh in the SCF "// &
                          "calculation, and the k mesh used for the NLO susceptibility "// &
                          "separately. Higher-order susceptibilities "// &
                          "usually require large optical k meshes.", &
                          n_keywords=10, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description= &
                          "Activates the nonlinear optical susceptibility calculation.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORDER", &
                          description= &
                          "Order n of the electric-field response $χ^{(n)}"//"$; "// &
                          "supported values are 2, 3, 4, and 5. Exactly "// &
                          "FREQ_1 through FREQ_n must be set.", &
                          usage="ORDER 2", type_of_var=integer_t, n_var=1, default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ELEMENTS", &
                          description= &
                          "Susceptibility elements to compute and print. Use -1 alone for "// &
                          "all elements (default), or a list of Cartesian index strings. "// &
                          "Each string has ORDER+1 letters X, Y, or Z: the output direction "// &
                          "followed by the input directions in FREQ_1 through FREQ_ORDER "// &
                          "order. For example, YXX selects $χ^{(2)}_{y;xx}"//"$, and XYYYX "// &
                          "selects $χ^{(4)}_{x;yyyx}"//"$. Uppercase and "// &
                          "lowercase are accepted."// &
                          " Repeat the keyword to continue a long list. Duplicate entries "// &
                          "are computed once. Only selected elements and their required "// &
                          "intermediate responses are evaluated; unselected "// &
                          "elements are omitted from the output. Do not "// &
                          "combine -1 with explicit elements.", &
                          usage="ELEMENTS YXX YYY", type_of_var=char_t, n_var=-1, &
                          repeats=.TRUE., default_c_vals=["-1"])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      DO i_frequency = 1, 5
         WRITE (frequency_name, '(A,I1)') "FREQ_", i_frequency
         CALL keyword_create(keyword, __LOCATION__, name=TRIM(frequency_name), &
                             description= &
                             "Input photon energy in eV. VALUE specifies one energy; RANGE "// &
                             "gives inclusive start, end, and positive step; SAME_AS and "// &
                             "MINUS_OF reference an earlier input frequency. "// &
                             "Independent ranges produce every combination. The sum "// &
                             "of input energies must be nonzero.", &
                             usage=TRIM(frequency_name)//" RANGE 0.5 1.0 0.01", &
                             type_of_var=char_t, n_var=-1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      END DO

      CALL keyword_create(keyword, __LOCATION__, name="T1_T2", &
                          description= &
                          "Common relaxation time $T_1=T_2=T$. Usually T₁ relaxes band "// &
                          "occupations (the diagonal density-matrix elements $ρ_{nn}"//"$) "// &
                          "towards equilibrium; T₂ describes decoherence, damping "// &
                          "off-diagonal elements $ρ_{mn}"//"$ to zero in the energy "// &
                          "eigenbasis. "// &
                          "Using one time gives $-(ρ-ρ_{eq})/T$ and removes the "// &
                          "basis-dependent distinction between populations and coherences "// &
                          "inside degenerate subspaces. This provides a simple, numerically "// &
                          "robust first approximation near band degeneracies and close "// &
                          "bands. It does not establish that the physical lifetimes are "// &
                          "equal. For relaxation-dependent singularities near the graphene "// &
                          "Dirac point, see Cheng, Vermeulen and Sipe (2015). The "// &
                          "broadening is $Γ=ℏ/T_2$. Must be positive; the default is 30 fs.", &
                          usage="T1_T2 [fs] 30", type_of_var=real_t, n_var=1, unit_str="fs", &
                          default_r_val=cp_unit_to_cp2k(value=30.0_dp, unit_str="fs"), &
                          citations=[Cheng2015])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="K_MESH", &
                          description= &
                          "Monkhorst-Pack-type uniform k mesh, shifted to include Γ, for the"// &
                          " Brillouin-zone integration required for $χ^{(n)}"//"$. Fractional "// &
                          "coordinates are k_i=j_i/N_i with j_i=0,...,N_i-1; no symmetry "// &
                          "reduction is used. This mesh can differ from the required k-point"// &
                          " SCF mesh. Higher-order nonlinear susceptibilities usually "// &
                          "require denser k meshes. Every entry must be "// &
                          "positive; use one in each nonperiodic direction.", &
                          usage="K_MESH 24 24 1", type_of_var=integer_t, n_var=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_nonlinear_optics_section

! **************************************************************************************************
!> \brief Read the activated input into nlo_env and reject invalid frequency combinations.
!> \param qs_env Converged QS environment containing the parsed input.
!> \param nlo_env Validated settings for the nonlinear response.
! **************************************************************************************************
   SUBROUTINE read_nlo_input(qs_env, nlo_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(nlo_type), INTENT(OUT)                        :: nlo_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_nlo_input'

      INTEGER                                            :: handle, i_case
      LOGICAL                                            :: zero_output
      REAL(KIND=dp)                                      :: energy_ev(5)
      TYPE(section_vals_type), POINTER                   :: input, section

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input)
      section => section_vals_get_subs_vals(input, 'PROPERTIES%NONLINEAR_OPTICS')
      CALL read_nlo_values(section, nlo_env)
      DO i_case = 1, nlo_env%n_frequency_combinations
         CALL nlo_frequency_case(nlo_env%frequency, nlo_env%order, i_case, energy_ev, zero_output)
         IF (zero_output) CPABORT('NLO output frequency must be nonzero.')
      END DO

      CALL timestop(handle)

   END SUBROUTINE read_nlo_input

! **************************************************************************************************
!> \brief Read and validate the NONLINEAR_OPTICS input before the first response calculation.
!> \param section Activated NONLINEAR_OPTICS input section.
!> \param nlo_env Validated response settings and selected elements.
! **************************************************************************************************
   SUBROUTINE read_nlo_values(section, nlo_env)
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(nlo_type), INTENT(OUT)                        :: nlo_env

      CHARACTER(len=6)                                   :: frequency_name
      CHARACTER(len=default_string_length), POINTER      :: words(:)
      INTEGER                                            :: i_frequency
      INTEGER, POINTER                                   :: mesh(:)
      TYPE(val_type), POINTER                            :: input_value

      nlo_env = nlo_type()
      IF (.NOT. ASSOCIATED(section)) CPABORT('NLO input section is missing.')
      NULLIFY (input_value, mesh, words)
      CALL section_vals_val_get(section, 'ORDER', i_val=nlo_env%order)
      IF (nlo_env%order < 2 .OR. nlo_env%order > 5) THEN
         CPABORT('NLO ORDER must be between 2 and 5.')
      END IF
      CALL read_nlo_elements(section, nlo_env)
      CALL section_vals_val_get(section, 'T1_T2', r_val=nlo_env%t1_t2_au)
      IF (.NOT. ieee_is_finite(nlo_env%t1_t2_au) .OR. nlo_env%t1_t2_au <= 0.0_dp) THEN
         CPABORT('NONLINEAR_OPTICS T1_T2 must be a positive, finite time.')
      END IF
      CALL section_vals_val_get(section, 'K_MESH', val=input_value)
      IF (.NOT. ASSOCIATED(input_value)) CPABORT('NONLINEAR_OPTICS requires K_MESH.')
      CALL section_vals_val_get(section, 'K_MESH', i_vals=mesh)
      IF (.NOT. ASSOCIATED(mesh)) CPABORT('NLO requires three K_MESH entries.')
      IF (SIZE(mesh) /= 3) CPABORT('NLO requires three K_MESH entries.')
      IF (ANY(mesh <= 0)) CPABORT('NONLINEAR_OPTICS K_MESH entries must be positive.')
      nlo_env%k_mesh = mesh
      DO i_frequency = 1, 5
         WRITE (frequency_name, '(A,I1)') 'FREQ_', i_frequency
         NULLIFY (input_value)
         CALL section_vals_val_get(section, frequency_name, val=input_value)
         IF (i_frequency > nlo_env%order) THEN
            IF (ASSOCIATED(input_value)) THEN
               CPABORT('Do not set input frequencies beyond NONLINEAR_OPTICS ORDER.')
            END IF
            CYCLE
         END IF
         IF (.NOT. ASSOCIATED(input_value)) CPABORT('Set all FREQ_1 through FREQ_ORDER.')
         NULLIFY (words)
         CALL section_vals_val_get(section, frequency_name, c_vals=words)
         IF (.NOT. ASSOCIATED(words)) CPABORT('NLO frequency specification is empty.')
         CALL parse_nlo_frequency(words, i_frequency, nlo_env%frequency(i_frequency))
      END DO
      ! Number of independent frequency tuples: product_j N(FREQ_j).
      CALL nlo_frequency_case_count(nlo_env%frequency, nlo_env%order, &
                                    nlo_env%n_frequency_combinations)

   END SUBROUTINE read_nlo_values

! **************************************************************************************************
!> \brief Combine repeated ELEMENTS lists without imposing a limit on their length.
!> \param section Input section containing the ELEMENTS keywords.
!> \param nlo_env Response settings receiving the selected elements.
! **************************************************************************************************
   SUBROUTINE read_nlo_elements(section, nlo_env)
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      CHARACTER(len=default_string_length), POINTER      :: words(:)
      INTEGER                                            :: i_rep, n_rep
      LOGICAL, ALLOCATABLE                               :: elements(:, :)

      ALLOCATE (nlo_env%elements(3, 3**nlo_env%order), elements(3, 3**nlo_env%order))
      nlo_env%elements = .FALSE.
      CALL section_vals_val_get(section, 'ELEMENTS', n_rep_val=n_rep)
      DO i_rep = 1, MAX(1, n_rep)
         NULLIFY (words)
         CALL section_vals_val_get(section, 'ELEMENTS', i_rep_val=i_rep, c_vals=words)
         IF (.NOT. ASSOCIATED(words)) CPABORT('NLO ELEMENTS list is empty.')
         IF (ANY(words == '-1') .AND. n_rep > 1) THEN
            CPABORT('NLO ELEMENTS: -1 cannot be repeated.')
         END IF
         CALL parse_nlo_elements(words, nlo_env%order, elements)
         nlo_env%elements(:, :) = nlo_env%elements .OR. elements
      END DO

   END SUBROUTINE read_nlo_elements

! **************************************************************************************************
!> \brief Decode Cartesian output/input strings into the selected tensor elements.
!> \param tokens Cartesian tensor-index strings from an ELEMENTS keyword.
!> \param order Number of incident electric fields.
!> \param elements Unique, validated tensor-index selections.
!> \note elements(c,1+sum_j(a_j-1)*3**(j-1)) selects χ_(c;a_1,...,a_order).
! **************************************************************************************************
   SUBROUTINE parse_nlo_elements(tokens, order, elements)
      CHARACTER(len=*), INTENT(IN)                       :: tokens(:)
      INTEGER, INTENT(IN)                                :: order
      LOGICAL, INTENT(OUT)                               :: elements(:, :)

      CHARACTER(len=6)                                   :: element
      INTEGER                                            :: i_code, i_direction, i_letter, i_out, &
                                                            i_token

      elements = .FALSE.
      IF (order < 2 .OR. order > 5) THEN
         CPABORT('NONLINEAR_OPTICS ELEMENTS requires ORDER 2 through 5.')
      END IF
      IF (ANY(SHAPE(elements) /= [3, 3**order]) .OR. SIZE(tokens) == 0) THEN
         CPABORT('NLO ELEMENTS list is empty.')
      END IF
      IF (ANY(tokens == '-1')) THEN
         IF (SIZE(tokens) /= 1) CPABORT('NLO ELEMENTS: -1 must appear alone.')
         elements = .TRUE.
         RETURN
      END IF
      DO i_token = 1, SIZE(tokens)
         IF (LEN_TRIM(tokens(i_token)) /= order + 1) THEN
            CPABORT('NLO ELEMENTS: use ORDER+1 letters, output first.')
         END IF
         element = tokens(i_token)
         CALL uppercase(element)
         i_out = 0
         i_code = 1
         DO i_letter = 1, order + 1
            i_direction = INDEX('XYZ', element(i_letter:i_letter))
            IF (i_direction == 0) THEN
               CPABORT('NONLINEAR_OPTICS ELEMENTS: use only X, Y, Z in each element.')
            END IF
            IF (i_letter == 1) THEN
               i_out = i_direction
            ELSE
               i_code = i_code + (i_direction - 1)*3**(i_letter - 2)
            END IF
         END DO
         elements(i_out, i_code) = .TRUE.
      END DO

   END SUBROUTINE parse_nlo_elements

! **************************************************************************************************
!> \brief Parse VALUE, RANGE, SAME_AS or MINUS_OF for input frequency i_frequency.
!> \param tokens Values of the incident-frequency keyword.
!> \param i_frequency One-based index of the incident frequency.
!> \param frequency Parsed single, ranged, or linked incident frequency.
! **************************************************************************************************
   SUBROUTINE parse_nlo_frequency(tokens, i_frequency, frequency)
      CHARACTER(len=*), INTENT(IN)                       :: tokens(:)
      INTEGER, INTENT(IN)                                :: i_frequency
      TYPE(nlo_frequency_type), INTENT(OUT)              :: frequency

      CHARACTER(len=32)                                  :: mode, reference_name
      INTEGER                                            :: i_reference, io_error
      REAL(KIND=dp)                                      :: end_ev, interval_count

      frequency = nlo_frequency_type()
      IF (i_frequency < 1 .OR. i_frequency > 5 .OR. SIZE(tokens) < 1) THEN
         CPABORT('Invalid FREQ index or empty frequency specification.')
      END IF
      mode = tokens(1)
      CALL uppercase(mode)
      SELECT CASE (TRIM(mode))
      CASE ('VALUE')
         IF (SIZE(tokens) /= 2) CPABORT('FREQ VALUE needs one energy in eV.')
         READ (tokens(2), *, IOSTAT=io_error) frequency%start_ev
         IF (io_error /= 0) CPABORT('FREQ VALUE requires a finite photon energy in eV.')
         IF (.NOT. ieee_is_finite(frequency%start_ev)) THEN
            CPABORT('FREQ VALUE requires a finite photon energy in eV.')
         END IF
         frequency%n_values = 1
      CASE ('RANGE')
         IF (SIZE(tokens) /= 4) CPABORT('FREQ RANGE requires start, end and step in eV.')
         READ (tokens(2), *, IOSTAT=io_error) frequency%start_ev
         IF (io_error /= 0) CPABORT('Invalid FREQ RANGE start energy.')
         READ (tokens(3), *, IOSTAT=io_error) end_ev
         IF (io_error /= 0) CPABORT('Invalid FREQ RANGE end energy.')
         READ (tokens(4), *, IOSTAT=io_error) frequency%step_ev
         IF (io_error /= 0) CPABORT('Invalid FREQ RANGE step energy.')
         IF (.NOT. ieee_is_finite(frequency%start_ev) .OR. .NOT. ieee_is_finite(end_ev) .OR. &
             .NOT. ieee_is_finite(frequency%step_ev)) THEN
            CPABORT('FREQ RANGE energies and step must be finite.')
         END IF
         IF (frequency%step_ev <= 0.0_dp .OR. end_ev < frequency%start_ev) THEN
            CPABORT('FREQ RANGE requires positive step and end >= start.')
         END IF
         interval_count = (end_ev - frequency%start_ev)/frequency%step_ev
         IF (.NOT. ieee_is_finite(interval_count) &
             .OR. interval_count >= REAL(HUGE(frequency%n_values) - 1, KIND=dp)) THEN
            CPABORT('FREQ RANGE contains too many values.')
         END IF
         frequency%n_values = 1 + FLOOR(interval_count + &
                                        32.0_dp*EPSILON(interval_count)*MAX(1.0_dp, interval_count))
      CASE ('SAME_AS', 'MINUS_OF')
         IF (SIZE(tokens) /= 2) CPABORT('FREQ relation needs one earlier index.')
         reference_name = tokens(2)
         CALL uppercase(reference_name)
         IF (LEN_TRIM(reference_name) /= 6 .OR. reference_name(1:5) /= 'FREQ_') THEN
            CPABORT('Frequency reference must be FREQ_1, FREQ_2, etc.')
         END IF
         READ (reference_name(6:6), *, IOSTAT=io_error) i_reference
         IF (io_error /= 0) CPABORT('Invalid frequency reference index.')
         IF (i_reference < 1 .OR. i_reference >= i_frequency) THEN
            CPABORT('Frequency references must point to an earlier FREQ.')
         END IF
         frequency%reference = i_reference
         IF (TRIM(mode) == 'MINUS_OF') frequency%sign = -1
      CASE DEFAULT
         CPABORT('FREQ must use VALUE, RANGE, SAME_AS or MINUS_OF.')
      END SELECT

   END SUBROUTINE parse_nlo_frequency

! **************************************************************************************************
!> \brief Count independent frequency combinations without materializing the Cartesian product.
!> \param frequencies Parsed incident-frequency definitions.
!> \param order Number of active incident frequencies.
!> \param n_frequency_combinations Product of independent range lengths.
! **************************************************************************************************
   SUBROUTINE nlo_frequency_case_count(frequencies, order, n_frequency_combinations)
      TYPE(nlo_frequency_type), INTENT(IN)               :: frequencies(5)
      INTEGER, INTENT(IN)                                :: order
      INTEGER, INTENT(OUT)                               :: n_frequency_combinations

      INTEGER                                            :: i_frequency

      n_frequency_combinations = 1
      IF (order < 2 .OR. order > 5) CPABORT('NLO ORDER must be between 2 and 5.')
      DO i_frequency = 1, order
         IF (frequencies(i_frequency)%reference > 0) THEN
            IF (frequencies(i_frequency)%reference >= i_frequency) THEN
               CPABORT('Frequency reference must point to an earlier FREQ.')
            END IF
         ELSE
            IF (frequencies(i_frequency)%n_values < 1) THEN
               CPABORT('Every input frequency through ORDER must be specified.')
            END IF
            IF (n_frequency_combinations > &
                HUGE(n_frequency_combinations)/frequencies(i_frequency)%n_values) THEN
               CPABORT('Too many independent frequency combinations.')
            END IF
            n_frequency_combinations = n_frequency_combinations*frequencies(i_frequency)%n_values
         END IF
      END DO

   END SUBROUTINE nlo_frequency_case_count

! **************************************************************************************************
!> \brief Resolve one Cartesian-product case, preserving SAME_AS and MINUS_OF links.
!> \param frequencies Parsed incident-frequency definitions.
!> \param order Number of active incident frequencies.
!> \param i_case One-based index of the requested range combination.
!> \param energy_ev Resolved incident photon energies in eV.
!> \param zero_output Whether the sum of incident photon energies is zero.
! **************************************************************************************************
   SUBROUTINE nlo_frequency_case(frequencies, order, i_case, energy_ev, zero_output)
      TYPE(nlo_frequency_type), INTENT(IN)               :: frequencies(5)
      INTEGER, INTENT(IN)                                :: order, i_case
      REAL(KIND=dp), INTENT(OUT)                         :: energy_ev(5)
      LOGICAL, INTENT(OUT)                               :: zero_output

      INTEGER                                            :: case_offset, i_frequency, i_value, &
                                                            n_frequency_combinations
      REAL(KIND=dp)                                      :: scale

      energy_ev = 0.0_dp
      zero_output = .FALSE.
      CALL nlo_frequency_case_count(frequencies, order, n_frequency_combinations)
      IF (i_case < 1 .OR. i_case > n_frequency_combinations) THEN
         CPABORT('Frequency case index is out of range.')
      END IF
      case_offset = i_case - 1
      DO i_frequency = 1, order
         IF (frequencies(i_frequency)%reference > 0) THEN
            energy_ev(i_frequency) &
               = frequencies(i_frequency)%sign*energy_ev(frequencies(i_frequency)%reference)
         ELSE
            i_value = MOD(case_offset, frequencies(i_frequency)%n_values)
            case_offset = case_offset/frequencies(i_frequency)%n_values
            energy_ev(i_frequency) = frequencies(i_frequency)%start_ev &
                                     + i_value*frequencies(i_frequency)%step_ev
         END IF
      END DO
      scale = MAX(1.0_dp, SUM(ABS(energy_ev(1:order))))
      zero_output = ABS(SUM(energy_ev(1:order))) <= 64.0_dp*EPSILON(scale)*scale

   END SUBROUTINE nlo_frequency_case

END MODULE nlo_read_input
