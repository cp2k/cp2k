!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief builds the input structure for optimize_basis
!> \par History
!>      03.2012 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************
MODULE input_optimize_basis

   USE cp_units, ONLY: cp_unit_to_cp2k
   USE input_constants, ONLY: do_opt_all, &
                              do_opt_coeff, &
                              do_opt_exps, &
                              do_opt_none
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type
   USE input_val_types, ONLY: char_t, &
                              integer_t, &
                              real_t
   USE kinds, ONLY: dp
   USE string_utilities, ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_optimize_basis'
   PUBLIC :: create_optimize_basis_section

CONTAINS

! **************************************************************************************************
!> \brief creates the optimize_basis section
!> \param section ...
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE create_optimize_basis_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPTIMIZE_BASIS", &
                          description="describes a basis optimization job, in which an ADMM like approach is used to"// &
                          " find the best exponents and/or coefficients to match a given training set.", &
                          repeats=.FALSE.)
      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="BASIS_TEMPLATE_FILE", &
                          description="Name of the basis set file, containing the structure of the new basis set", &
                          usage="BASIS_TEMPLATE_FILE <FILENAME>", &
                          type_of_var=char_t, repeats=.FALSE., &
                          default_c_val="BASIS_SET", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BASIS_WORK_FILE", &
                          description="Name of the basis set file which is created to be read as initial guess", &
                          usage="BASIS_WORK_FILE <FILENAME>", &
                          type_of_var=char_t, repeats=.FALSE., &
                          default_c_val="BASIS_WORK_FILE", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BASIS_OUTPUT_FILE", &
                          description="Name of the basis set file containing the optimized basis", &
                          usage="BASIS_OUTPUT_FILE <FILENAME>", &
                          type_of_var=char_t, repeats=.FALSE., &
                          default_c_val="BASIS_OUTPUT_FILE", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WRITE_FREQUENCY", &
                          description="Frequency at which the intermediate results should be written", &
                          usage="WRITE_FREQUENCY 1000", &
                          default_i_val=5000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_CONDITION_NUMBER", &
                          description="Determines whether condition number should be part of optimization or not", &
                          usage="USE_CONDITION_NUMBER", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="BASIS_COMBINATIONS", &
         description="If multiple atomic kinds are fitted at the same time, this keyword "// &
         "allows to specify which basis sets should be used together in optimization (underived set ID=0). "// &
         "If skipped all combinations are used. The order is taken as the kinds and sets are specified in the input", &
         repeats=.TRUE., &
         usage="BASIS_COMBINATIONS SET_ID(KIND1) SET_ID(KIND2) ... ", type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="RESIDUUM_WEIGHT", &
         description="This keyword allows to give different weight factors to the "// &
         "residuum of the different basis combinations. "// &
         "The first entry corresponds to the original basis sets. Every further value is assigned to the combinations "// &
         "in the order given for BASIS_COMBINATIONS.", &
         repeats=.TRUE., &
         usage="RESIDUUM_WEIGHT REAL ", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="CONDITION_WEIGHT", &
         description="This keyword allows to give different weight factors to the "// &
         "condition number of different basis combinations (LOG(cond) is used). "// &
         "The first entry corresponds to the original basis sets. Every further value is assigned to the combinations "// &
         "in the order given for BASIS_COMBINATIONS.", &
         repeats=.TRUE., &
         usage="CONDITION_WEIGHT REAL ", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GROUP_PARTITION", &
                          description="Allows the specification of the group mpi group sizes in parallel "// &
                          "runs. If less Groups than tasks are speciefied, consecutive calculations "// &
                          "Will be assigned to one group (derived basis sets and then training sets) "// &
                          "If keyword is skipped, equal group sizes will be generated trying to fit all calculations.", &
                          repeats=.TRUE., &
                          usage="GROUP_PARTITION INT INT ... ", type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_fit_kinds_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_frontier_orbitals_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_frontier_orbital_screening_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_training_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_powell_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_optimize_basis_section

! **************************************************************************************************
!> \brief Creates the optional frontier-orbital objective section.
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_frontier_orbitals_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=2), PARAMETER                        :: paragraph_break = NEW_LINE("A")//NEW_LINE("A")
      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create( &
         section, __LOCATION__, name="FRONTIER_ORBITALS", &
         description="Optimizes a compact basis set to reproduce the occupied and low-energy virtual "// &
         "Kohn-Sham states obtained from a very large reference basis set. For each training system, "// &
         "a quantum-mechanical DFT reference calculation using this basis must be provided through "// &
         "TRAINING_FILES, including a converged reference wavefunction. The Kohn-Sham Hamiltonian from "// &
         "the reference calculation is frozen and "// &
         "projected into every trial compact basis. Consequently, all candidate orbital subspaces and "// &
         "energies, including the HOMO-LUMO gap, are evaluated relative to the reference Hamiltonian. "// &
         "When optimizing the compact basis against the large reference basis, the following loss "// &
         "function is minimized:"//paragraph_break// &
         "$$\mathcal{L}=\alpha_{\mathrm{occ}}\mathcal{L}_{\mathrm{occ}}+ "// &
         "\alpha_{\mathrm{vir}}\mathcal{L}_{\mathrm{vir}}+ "// &
         "\alpha_{\mathrm{empty}}\mathcal{L}_{\mathrm{empty}}+ "// &
         "\alpha_{\mathrm{gap}}\mathcal{L}_{\mathrm{gap}}+ "// &
         "\alpha_{\mathrm{coeff}}\mathcal{L}_{\mathrm{coeff}}+ "// &
         "\alpha_\kappa\mathcal{L}_\kappa.$$"//paragraph_break// &
         "Below, all loss functions $\mathcal{L}$ are defined. $R$ and $C$ label reference and "// &
         "candidate quantities, $i$ and $j$ label occupied orbitals, and $a$ and $b$ label virtual "// &
         "orbitals. Their overlaps are $M_{ij}=\langle\psi_i^R\vert\psi_j^C\rangle$ and "// &
         "$M_{ab}=\langle\psi_a^R\vert\psi_b^C\rangle$, respectively. "// &
         "$N_{\mathrm{occ}}$ is the number of occupied spatial orbitals. "// &
         "The occupied-subspace loss is"//paragraph_break// &
         "$$\mathcal{L}_{\mathrm{occ}}=1-\frac{1}{N_{\mathrm{occ}}} "// &
         "\sum_{i,j\in\mathrm{occ}}\lvert M_{ij}\rvert^2.$$"//paragraph_break// &
         "It penalizes a poor representation of the reference occupied orbitals by the occupied "// &
         "orbitals of the compact candidate basis. "// &
         "Low-energy virtual orbitals are selected smoothly using"//paragraph_break// &
         "$$w_a^X=\left[1+\exp\left( "// &
         "\frac{\epsilon_a^X-\epsilon_{\mathrm{LUMO}}^X-E_{\mathrm{cut}}}{\Delta E} "// &
         "\right)\right]^{-1},\qquad X\in\{R,C\}.$$"//paragraph_break// &
         "Define $N_X=\sum_{a\in\mathrm{vir}}(w_a^X)^2$. The virtual-subspace loss can be written "// &
         "directly as the normalized squared distance between the two energy-weighted virtual-space "// &
         "operators:"// &
         paragraph_break// &
         "$$\begin{aligned} "// &
         "\mathcal{L}_{\mathrm{vir}}&=\frac{1}{2N_R}\left\| "// &
         "\sum_{a\in\mathrm{vir}}w_a^R\lvert\psi_a^R\rangle\langle\psi_a^R\rvert- "// &
         "\sum_{b\in\mathrm{vir}}w_b^C\lvert\psi_b^C\rangle\langle\psi_b^C\rvert "// &
         "\right\|_F^2\\ "// &
         "&=\frac{N_R+N_C-2\sum_{a,b\in\mathrm{vir}} "// &
         "w_a^R w_b^C\lvert M_{ab}\rvert^2}{2N_R}. "// &
         "\end{aligned}$$"//paragraph_break// &
         "The first operator represents the reference low-energy virtual space and the second represents "// &
         "the candidate low-energy virtual space. Here, $\lVert\cdot\rVert_F^2$ is the squared Frobenius "// &
         "norm, which measures their squared matrix distance. Thus, $\mathcal{L}_{\mathrm{vir}}$ penalizes "// &
         "differences between the low-energy virtual orbitals of the compact candidate basis and the "// &
         "reference space; it is zero when these spaces coincide. The compact candidate basis can contain "// &
         "empty orbitals outside the "// &
         "selected low-energy interval. To test whether the reference low-energy virtual orbitals can be "// &
         "represented in the entire candidate empty space, the empty-subspace loss is"//paragraph_break// &
         "$$\mathcal{L}_{\mathrm{empty}}=1- "// &
         "\frac{\sum_{a,b\in\mathrm{vir}}(w_a^R)^2\lvert M_{ab}\rvert^2}{N_R}. "// &
         "$$"//paragraph_break// &
         "It penalizes reference low-energy virtual orbitals that cannot be represented by any empty "// &
         "orbital of the compact candidate basis. "// &
         "With $E_{\mathrm{gap}}^X=\epsilon_{\mathrm{LUMO}}^X-\epsilon_{\mathrm{HOMO}}^X$, "// &
         "the gap loss is"//paragraph_break// &
         "$$\mathcal{L}_{\mathrm{gap}}=\left( "// &
         "\frac{E_{\mathrm{gap}}^C-E_{\mathrm{gap}}^R}{E_{\mathrm{scale}}} "// &
         "\right)^2.$$"//paragraph_break// &
         "It penalizes deviations of the candidate HOMO-LUMO gap from the reference gap."// &
         paragraph_break// &
         "During the optimization, the electron density and the corresponding Kohn-Sham Hamiltonian "// &
         "are held fixed at the values obtained with the large reference basis set. Therefore, the "// &
         "candidate gap reported in the optimization output is the fixed-density gap of the reference "// &
         "Hamiltonian projected into the compact candidate basis. A separate self-consistent Kohn-Sham "// &
         "DFT calculation with the optimized compact basis can produce a different electron density and, "// &
         "consequently, a different HOMO-LUMO gap."//paragraph_break// &
         "For the "// &
         "optimized and initial contraction coefficients of the compact basis set, $c$ and $c_0$, "// &
         "the coefficient loss is"//paragraph_break// &
         "$$\mathcal{L}_{\mathrm{coeff}}= "// &
         "\frac{\sum_p(c_p-c_{0,p})^2}{\sum_p c_{0,p}^2}.$$"//paragraph_break// &
         "It penalizes large deviations from the coefficients of the initially supplied compact "// &
         "basis set. Finally, the conditioning loss is"//paragraph_break// &
         "$$\mathcal{L}_\kappa=\log_{10}\kappa(S_C),\qquad "// &
         "\kappa(S_C)=\frac{\lambda_{\max}(S_C)}{\lambda_{\min}(S_C)}.$$"//paragraph_break// &
         "It penalizes an ill-conditioned candidate overlap matrix. The weights "// &
         "$\alpha_{\mathrm{occ}}$, $\alpha_{\mathrm{vir}}$, $\alpha_{\mathrm{empty}}$, "// &
         "$\alpha_{\mathrm{gap}}$, and $\alpha_{\mathrm{coeff}}$ are set below; "// &
         "CONDITION_WEIGHT supplies $\alpha_\kappa$."//paragraph_break// &
         "The loss weights $\alpha_{\mathrm{occ}}$, $\alpha_{\mathrm{vir}}$, "// &
         "$\alpha_{\mathrm{empty}}$, $\alpha_{\mathrm{gap}}$, "// &
         "$\alpha_{\mathrm{coeff}}$, and $\alpha_\kappa$ can have a strong and "// &
         "system-dependent effect on the accuracy of the optimized basis set. In particular, the user "// &
         "has to test the new basis set in a self-consistent DFT calculation, which also uses the "// &
         "self-consistent electron density computed with the optimized basis set (during the "// &
         "optimization, the electron density is kept fixed at the density obtained with the large "// &
         "reference basis set). Several combinations of loss weights should be tested, and the optimal "// &
         "basis should be selected based on validation with a self-consistent DFT calculation using the "// &
         "new basis set (for example, by focusing on the resulting DFT gap, GW gap, or GW-BSE excitation "// &
         "energy). Useful candidate parameter combinations include:"//paragraph_break// &
         "1. OCCUPIED_WEIGHT 1, VIRTUAL_WEIGHT 1, EMPTY_OVERLAP_WEIGHT 1, GAP_WEIGHT 10, "// &
         "GAP_ENERGY_SCALE [eV] 1, CONDITION_WEIGHT 0.1, COEFFICIENT_WEIGHT 0.001"//NEW_LINE("A")// &
         "2. OCCUPIED_WEIGHT 10, VIRTUAL_WEIGHT 10, EMPTY_OVERLAP_WEIGHT 10, GAP_WEIGHT 10, "// &
         "GAP_ENERGY_SCALE [eV] 1, CONDITION_WEIGHT 0.001, COEFFICIENT_WEIGHT 0.001"//NEW_LINE("A")// &
         "3. OCCUPIED_WEIGHT 30, VIRTUAL_WEIGHT 30, EMPTY_OVERLAP_WEIGHT 30, GAP_WEIGHT 10, "// &
         "GAP_ENERGY_SCALE [eV] 1, CONDITION_WEIGHT 0.001, COEFFICIENT_WEIGHT 0.001"//NEW_LINE("A")// &
         "4. OCCUPIED_WEIGHT 10, VIRTUAL_WEIGHT 30, EMPTY_OVERLAP_WEIGHT 30, GAP_WEIGHT 10, "// &
         "GAP_ENERGY_SCALE [eV] 1, CONDITION_WEIGHT 0.001, COEFFICIENT_WEIGHT 0.01", &
         repeats=.FALSE.)

      CALL keyword_create( &
         keyword, __LOCATION__, name="OCCUPIED_WEIGHT", &
         description="Prefactor $\alpha_{\mathrm{occ}}$ multiplying the occupied-subspace loss "// &
         "$\mathcal{L}_{\mathrm{occ}}$.", &
         usage="OCCUPIED_WEIGHT 1.0", type_of_var=real_t, default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="VIRTUAL_WEIGHT", &
         description="Prefactor $\alpha_{\mathrm{vir}}$ multiplying the virtual-subspace loss "// &
         "$\mathcal{L}_{\mathrm{vir}}$.", &
         usage="VIRTUAL_WEIGHT 1.0", type_of_var=real_t, default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EMPTY_OVERLAP_WEIGHT", &
         description="Prefactor $\alpha_{\mathrm{empty}}$ multiplying "// &
         "$\mathcal{L}_{\mathrm{empty}}$, the loss of reference low-energy virtual-state "// &
         "completeness in the candidate virtual subspace.", &
         usage="EMPTY_OVERLAP_WEIGHT 1.0", type_of_var=real_t, default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="GAP_WEIGHT", &
         description="Prefactor $\alpha_{\mathrm{gap}}$ multiplying the normalized HOMO-LUMO gap loss "// &
         "$\mathcal{L}_{\mathrm{gap}}$. Because $\mathcal{L}_{\mathrm{gap}}$ is proportional to "// &
         "$1/E_{\mathrm{scale}}^2$, reducing GAP_ENERGY_SCALE by a factor "// &
         "of 100 requires reducing GAP_WEIGHT by a factor of 10000 to retain the same balance.", &
         usage="GAP_WEIGHT 1.0E-3", type_of_var=real_t, default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="COEFFICIENT_WEIGHT", &
         description="Prefactor $\alpha_{\mathrm{coeff}}$ multiplying the normalized squared change "// &
         "$\mathcal{L}_{\mathrm{coeff}}$ of the optimized contraction coefficients from their "// &
         "initial values.", &
         usage="COEFFICIENT_WEIGHT 1.0E-3", type_of_var=real_t, default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="GAP_ENERGY_SCALE", &
         description="$E_{\mathrm{scale}}$ in the gap loss. A gap error of this magnitude gives "// &
         "$\mathcal{L}_{\mathrm{gap}}=1$. The default is 0.01 eV.", &
         usage="GAP_ENERGY_SCALE [eV] 0.01", type_of_var=real_t, unit_str="eV", &
         default_r_val=cp_unit_to_cp2k(value=0.01_dp, unit_str="eV"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="VIRTUAL_ENERGY_CUTOFF", &
         description="$E_{\mathrm{cut}}$, the virtual-orbital energy window measured "// &
         "from the corresponding LUMO. The default is 3 eV.", &
         usage="VIRTUAL_ENERGY_CUTOFF [eV] 3.0", type_of_var=real_t, unit_str="eV", &
         default_r_val=cp_unit_to_cp2k(value=3.0_dp, unit_str="eV"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="VIRTUAL_ENERGY_SMOOTHING", &
         description="$\Delta E$, the smoothing width at the outer boundary of the "// &
         "virtual-orbital energy window. It must be positive; the default is 0.2 eV.", &
         usage="VIRTUAL_ENERGY_SMOOTHING [eV] 0.2", type_of_var=real_t, unit_str="eV", &
         default_r_val=cp_unit_to_cp2k(value=0.2_dp, unit_str="eV"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_frontier_orbitals_section

! **************************************************************************************************
!> \brief Creates the optional serial frontier-orbital screening section.
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_frontier_orbital_screening_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create( &
         section, __LOCATION__, name="FRONTIER_ORBITAL_SCREENING", &
         description="Optimize a basis set to match frontier orbitals of a reference basis "// &
         "(frontier orbitals: orbitals close to the HOMO and close to the LUMO). This section uses the "// &
         "optimization described in OPTIMIZE_BASIS%FRONTIER_ORBITALS, where the complete optimization "// &
         "theory is given. The optimization contains several parameters, in particular the prefactors "// &
         "of its different loss functions, and the resulting optimized basis set can depend strongly "// &
         "on these prefactors. This section runs up to five pre-tabulated sets of optimization "// &
         "parameters. All optimized basis sets are reported together with their self-consistent "// &
         "HOMO-LUMO gaps and total energies. The optimization uses the fixed electron density obtained "// &
         "with a large reference basis, and CP2K also converges an SCF calculation with every optimized "// &
         "basis set and reports the resulting HOMO-LUMO gap and total energy. No training calculation "// &
         "or previous DFT calculation needs to be read. All DFT calculations are performed automatically "// &
         "inside CP2K when this section is activated, using the DFT input parameters supplied in the "// &
         "FORCE_EVAL section of the input file. All optimized basis sets are written to "// &
         "BASIS_OUTPUT_FILE with suffixes such as _1 and _2.", &
         repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="NUMBER_OF_OPTIMIZATIONS", &
                          description="Number of pre-tabulated sets of optimization parameters.", &
                          usage="NUMBER_OF_OPTIMIZATIONS 5", &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_frontier_orbital_screening_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_fit_kinds_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword, subsection)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="FIT_KIND", &
                          description="specicifies the atomic kinds to be fitted and the basis"// &
                          " sets associated with the kind.", &
                          repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="The name of the kind described in this section.", &
                          usage="H", default_c_val="DEFAULT")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BASIS_SET", &
                          description="The name of the basis set for the kind. Has to be specified in BASIS_TEMPLATE_FILE.", &
                          usage="BASIS_SET H", default_c_val="DEFAULT")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INITIAL_DEGREES_OF_FREEDOM", &
                          description="Specifies the initial degrees of freedom in the basis optimization. "// &
                          "This can be used to make further specifications easier", &
                          usage="INITIAL_DEGREES_OF_FREEDOM ALL", &
                          enum_c_vals=s2a("ALL", "NONE", "COEFFICIENTS", "EXPONENTS"), &
                          enum_desc=s2a("Set all parameters in the basis to be variable.", &
                                        "Set all parameters in the basis to be fixed.", &
                                        "Set all coefficients in the basis set to be variable.", &
                                        "Set all exponents in the basis to be variable."), &
                          enum_i_vals=[do_opt_all, do_opt_none, do_opt_coeff, do_opt_exps], &
                          default_i_val=do_opt_coeff)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SWITCH_COEFF_STATE", &
                          description="Allows to switch the state of a given coefficient from current state "// &
                          "(varibale/fixed)) to the opposite state. The three integers indicate "// &
                          "the set number, the angular momentum i'th contraction and i'th coefficient", repeats=.TRUE., &
                          usage="SWITCH_COEFF_STATE SET L CONTRACTION IPGF", type_of_var=integer_t, n_var=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SWITCH_CONTRACTION_STATE", &
                          description="Allows to switch the state of a given contraction from current state "// &
                          "(varibale/fixed)) to the opposite state. The three integers indicate "// &
                          "the set number, the angular momentum and i'th contraction ", repeats=.TRUE., &
                          usage="SWITCH_CONTRACTION_STATE SET L CONTRACTION ", type_of_var=integer_t, n_var=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SWITCH_EXP_STATE", &
                          description="Allows to switch the state of a given exponent from current state "// &
                          "(varibale/fixed)) to the opposite state. The two integers indicate "// &
                          "the set number and i'th exponent", repeats=.TRUE., &
                          usage="SWITCH_EXP_STATE SET IEXP", type_of_var=integer_t, n_var=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SWITCH_SET_STATE", &
                          description="Allows to switch the states of in a set from current state "// &
                          "(varibale/fixed)) to the opposite state. The two integers indicate "// &
                          "the affected part (0=ALL,1=EXPS,2=COEFF) and i'th set", repeats=.TRUE., &
                          usage="SWITCH_SET_STATE SET IEXP", type_of_var=integer_t, n_var=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_constrain_exp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_derived_sets_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_fit_kinds_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_derived_sets_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DERIVED_BASIS_SETS", &
                          description="This section can be used to create subsets of a basis"// &
                          " which will be fitted at the same time. This is especially useful if connected"// &
                          " bsis sets e.g. TZVP, DZVP, SZV should be fitted.", &
                          repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="BASIS_SET_NAME", &
                          description="Defines the name of the derived basis set, which will be "// &
                          "automatically generated otherwise.", &
                          usage="BASIS_SET_NAME {word}", &
                          type_of_var=char_t, &
                          repeats=.FALSE., &
                          default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_SET", &
                          description="Specifies the reference basis ID which is used as template to create the new set. "// &
                          "The original basis has ID 0. All following sets are counted in order as specified in the Input."// &
                          " The descriptors always assume the structure of the input basis set.", &
                          repeats=.FALSE., usage="REFERENCE_SET INTEGER", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REMOVE_CONTRACTION", &
                          description="Can be used to remove a contraction from the reference basis set. "// &
                          "The contraction is speciefied by set number, angular momentum and number of contraction."// &
                          " The descriptors always assume the structure of the input basis set.", &
                          repeats=.TRUE., usage="REMOVE_CONTRACTION SET L ICONTRACTION", type_of_var=integer_t, n_var=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REMOVE_SET", &
                          description="Can be used to remove a set from the reference basis set. ", &
                          repeats=.TRUE., usage="REMOVE_SET SET", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_derived_sets_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_constrain_exp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CONSTRAIN_EXPONENTS", &
                          description="specicifies constraints for the exponents to be fitted."// &
                          " Only a single constraint can be applied to an exponent", &
                          repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="USE_EXP", &
                          description="Defines the exponent to be constraint. The two integers indicate "// &
                          "the set number and i'th exponent. The value -1 can be used to mark all sets/exponents in a set.", &
                          repeats=.FALSE., usage="USE_EXP SET IEXP", type_of_var=integer_t, n_var=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BOUNDARIES", &
                          description="Defines the boundaries to which the optimization is restricted."// &
                          " First value is the lower bound, second value is the upper bound.", &
                          repeats=.FALSE., usage="BOUNDARIES LOWER UPPER", type_of_var=real_t, n_var=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_VAR_FRACTION", &
                          description="Defines the maximum fractionr by which the exponent is allowed to vary."// &
                          " e.g. 0.5 allows the exp to vary by 0.5*exp in both directions.", &
                          repeats=.FALSE., usage="MAX_VAR_FRACTION REAL", type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_constrain_exp_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_training_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TRAINING_FILES", &
                          description="specicifies the location in which the files necessary for"// &
                          " fitting procedure are located. Each Training set needs a repetition of this section.", &
                          repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="DIRECTORY", &
                          description="the directory in which the files are placed", &
                          usage="DIRECTORY /my/path", &
                          default_lc_val=".")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INPUT_FILE_NAME", &
                          description="the filename of the input file used to run the original calculation", &
                          usage="INPUT_FILE_NAME my_input.inp", &
                          default_lc_val="input.inp")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_training_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_powell_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPTIMIZATION", &
                          description="sets the parameters for optimizition, output frequency and restarts", &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="ACCURACY", &
                          description="Final accuracy requested in optimization (RHOEND)", &
                          usage="ACCURACY 0.00001", &
                          default_r_val=1.e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STEP_SIZE", &
                          description="Initial step size for search algorithm (RHOBEG)", &
                          usage="STEP_SIZE 0.005", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_FUN", &
                          description="Maximum number of function evaluations", &
                          usage="MAX_FUN 1000", &
                          default_i_val=5000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_powell_section

END MODULE input_optimize_basis
