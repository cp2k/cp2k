!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input of the real-space resolution of the identity (RI-RS), shared by GW in
!>        PROPERTIES%BANDSTRUCTURE%GW and by RPA in DFT%XC%WF_CORRELATION%RI_RS_RPA
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE input_cp2k_ri_rs
!&<
   USE bibliography,                    ONLY: Duchemin2019, &
                                              Duchemin2021
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: ri_rs_d_lP_analytic, &
                                              ri_rs_d_lP_auto, &
                                              ri_rs_d_lP_numeric
   USE input_keyword_types,             ONLY: keyword_create, &
                                              keyword_release, &
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword, &
                                              section_add_subsection, &
                                              section_create, &
                                              section_release, &
                                              section_type
   USE input_val_types,                 ONLY: real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
!&>
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_ri_rs'

   PUBLIC :: add_eps_filter_keyword, &
             add_minimax_keywords, &
             add_ri_keywords, &
             create_auto_ri_section, &
             create_ri_rs_rpa_section, &
             create_ri_rs_section

CONTAINS

! **************************************************************************************************
!> \brief Adds the keywords of the minimax imaginary-time and imaginary-frequency grids
!> \param section GW or RI_RS_RPA input section
! **************************************************************************************************
   SUBROUTINE add_minimax_keywords(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NUM_TIME_FREQ_POINTS", &
                          description="Number of discrete points for the imaginary-time "// &
                          "grid and the imaginary-frequency grid. The more points, the more "// &
                          "precise is the calculation. Typically, 10 points are good "// &
                          "for 0.1 eV precision of band structures and molecular energy "// &
                          "levels, 20 points for 0.03 eV precision, "// &
                          "and 30 points for 0.01 eV precision, see Table I in "// &
                          "\[Wilhelm2021\]. The computation time of RPA and GW "// &
                          "increases linearly with `NUM_TIME_FREQ_POINTS`.", &
                          usage="NUM_TIME_FREQ_POINTS 30", &
                          default_i_val=30)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REGULARIZATION_MINIMAX", &
                          description="Parameter to regularize the Fourier transformation with minimax grids. "// &
                          "In case the parameter 0.0 is chosen, no regularization is performed.", &
                          usage="REGULARIZATION_MINIMAX 1.0E-4", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE add_minimax_keywords

! **************************************************************************************************
!> \brief Adds the keyword EPS_FILTER, the threshold of sparse matrix operations
!> \param section GW or RI_RS_RPA input section
! **************************************************************************************************
   SUBROUTINE add_eps_filter_keyword(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold for sparse operations on spatially local "// &
                          "matrices, including DBCSR matrix multiplication. EPS_FILTER "// &
                          "usually controls the accuracy and runtime of low-scaling RPA and GW: "// &
                          "a smaller value increases numerical precision and computational "// &
                          "cost.", &
                          usage="EPS_FILTER 1.0E-6", &
                          default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE add_eps_filter_keyword

! **************************************************************************************************
!> \brief Adds the keywords of the resolution of the identity: regularization and RI metric
!> \param section GW or RI_RS_RPA input section
! **************************************************************************************************
   SUBROUTINE add_ri_keywords(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REGULARIZATION_RI", &
                          description="Regularization applied to the RI fitting matrix $M_{PQ}$ of RPA "// &
                          "and GW. A negative value selects the context-dependent default: "// &
                          "zero for molecules, 1.0E-4 for periodic RI-RS RPA and GW and 1.0E-2 "// &
                          "for other periodic systems.", &
                          usage="REGULARIZATION_RI 1.0E-4", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS_RI", &
                          description="The cutoff radius (in Angstrom) for the truncated  "// &
                          "Coulomb operator. The larger the cutoff radius, the faster "// &
                          "converges the resolution of the identity (RI) with respect to the "// &
                          "RI basis set size. Larger cutoff radius means higher computational "// &
                          "cost. Not used in periodic RI-RS RPA and GW with k-points, where "// &
                          "RI_CUTOFF_ATTENUATION sets the RI metric.", &
                          usage="CUTOFF_RADIUS_RI 3.0", &
                          default_r_val=cp_unit_to_cp2k(value=3.0_dp, unit_str="angstrom"), &
                          type_of_var=real_t, unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RI_CUTOFF_ATTENUATION", &
                          description="Positive attenuation omega (in inverse bohr) of the "// &
                          "erfc(omega*r)/r RI metric used by periodic RI-RS RPA and GW with k-points. "// &
                          "The default is 0.3 bohr^-1. Other GW implementations use "// &
                          "the truncated Coulomb metric set by CUTOFF_RADIUS_RI.", &
                          usage="RI_CUTOFF_ATTENUATION 0.3", default_r_val=0.3_dp, &
                          type_of_var=real_t, unit_str="bohr^-1")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE add_ri_keywords

! **************************************************************************************************
!> \brief Creates the automatic RI basis set optimization section.
!> \param section AUTO_RI input section
! **************************************************************************************************
   SUBROUTINE create_auto_ri_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="AUTO_RI", &
                          description="Optimizes an RI basis for RI-RS RPA and GW, for molecules "// &
                          "and for periodic systems with k-points. The fitting target combines "// &
                          "occupied-empty transitions and exchange products.", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="RI_AO_RATIO", &
                          description="Ratio between the numbers of optimized RI and AO basis "// &
                          "functions. For example, RI_AO_RATIO 3 requires an RI basis with "// &
                          "three times as many functions as the AO basis. Compare values such "// &
                          "as 1, 1.5, 2, 3 and 4 until the RPA or GW energies are converged. "// &
                          "Larger values increase time and memory, but also increase the "// &
                          "numerical accuracy.", &
                          default_r_val=3.0_dp, usage="RI_AO_RATIO 3")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OCC_EMPTY_FRONTIER_ORBITAL_WINDOW", &
                          description="The exchange contribution to the RI optimization "// &
                          "objective includes GW levels within the energy window "// &
                          "$[\epsilon_{\mathrm{HOMO}}-w,\epsilon_{\mathrm{LUMO}}+w]$. "// &
                          "The default of 5 eV is usually sufficient, for example for "// &
                          "GW-BSE excitations up to 10 eV. GW only: RI_RS_RPA uses the default.", &
                          default_r_val=cp_unit_to_cp2k(value=5.0_dp, unit_str="eV"), &
                          unit_str="eV", usage="OCC_EMPTY_FRONTIER_ORBITAL_WINDOW [eV] 5.0")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_auto_ri_section

! **************************************************************************************************
!> \brief creates the input section collecting all parameters of the real-space resolution of
!>        the identity (RI-RS) used in RPA and GW
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_ri_rs_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_RS", &
                          description="Real-Space Resolution of Identity (RI-RS) method. This "// &
                          "approximation replaces the conventional 3-center RI integrals "// &
                          "$(\mu\nu|P)$ by a factorized representation on an atom-centered "// &
                          "real-space grid with points $\mathbf{r}_\ell$: "// &
                          "$(\mu\nu|P) \approx \sum_\ell \varphi_\mu(\mathbf{r}_\ell) "// &
                          "\varphi_\nu(\mathbf{r}_\ell) Z_{\ell P}$. "// &
                          "The coefficients $Z_{\ell P}$ combine the numerical integration "// &
                          "weights and the Coulomb potential of the auxiliary basis function "// &
                          "$P$ evaluated at grid point $\mathbf{r}_\ell$. To reduce the "// &
                          "computational cost, only grid points within the sphere $B^P$ around "// &
                          "atom $P$ are included, i.e. the points with "// &
                          "$|\mathbf{r}_\ell - \mathbf{R}_P| < R_c + r_P$. "// &
                          "Here, $r_P$ is the effective radius of the most diffuse RI "// &
                          "auxiliary Gaussian on atom $P$, at which the basis function "// &
                          "magnitude falls below a threshold $\delta$ "// &
                          "(currently controlled through EPS_FILTER). "// &
                          "This locality approximation yields a sparse representation of the "// &
                          "3-center integrals and enables reduced computational cost. "// &
                          "See details in \[Duchemin2019\].", &
                          citations=[Duchemin2019, Duchemin2021], &
                          n_keywords=12, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of RI-RS in the RPA or GW calculation. "// &
                          "An empty section, &RI_RS ... &END RI_RS, switches RI-RS on with "// &
                          "default parameters.", &
                          usage="&RI_RS .TRUE.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TIKHONOV", &
                          description="Regularization parameter $\alpha$ used to stabilize "// &
                          "the inversion of the grid-overlap matrix "// &
                          "$D$ in the RI-RS method. "// &
                          "See Equation (9) in \[Duchemin2019\]. A negative value selects "// &
                          "the context-dependent default: 1.0E-4 for periodic RI-RS RPA and GW with k-points "// &
                          "and 1.0E-8 otherwise.", &
                          usage="TIKHONOV 1.0E-8", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID_SELECT", &
                          description="Selection of the atom-centered grid type used "// &
                          "in RI-RS optimized by Duchemin and Blase. "// &
                          "(1) def2-TZVPP: Available for elements up to the fourth row "// &
                          "of the periodic table (see \[Duchemin2021\]). "// &
                          "(2) cc-pVTZ: Available for H, C, N, and O atoms "// &
                          "(see \[Duchemin2019\]). "// &
                          "(3) User-provided grids: per-element grid files supplied by the "// &
                          "user, read as ri_rs_grid/<Element><suffix> in the same format as "// &
                          "the built-in sets, with the suffix _rirs.ion.", &
                          usage="GRID_SELECT 1", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D_LP_COMPUTATION", &
                          description="Computation of "// &
                          "$d_{\ell P} = \sum_{\mu\nu} \varphi_\mu(\mathbf{r}_\ell) "// &
                          "\varphi_\nu(\mathbf{r}_\ell) (\mu\nu|P)$ in periodic RI-RS with "// &
                          "k-point sampling.", &
                          usage="D_LP_COMPUTATION AUTO", &
                          enum_c_vals=s2a("AUTO", "ANALYTIC", "NUMERIC"), &
                          enum_i_vals=[ri_rs_d_lP_auto, ri_rs_d_lP_analytic, ri_rs_d_lP_numeric], &
                          enum_desc=s2a("Estimates the computation time of ANALYTIC and NUMERIC "// &
                                        "and uses the faster one.", &
                                        "Analytic three-center integrals $(\mu\nu|P)$ with "// &
                                        "lattice sums; integrals related by a space-group "// &
                                        "operation are computed once.", &
                                        "Numerical integration on a uniform real-space grid. "// &
                                        "Not feasible for all-electron basis sets."), &
                          default_i_val=ri_rs_d_lP_auto)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS_RL_AO", &
                          description="Real-space cutoff radius (in Angstrom) for evaluating "// &
                          "the AO basis functions on the RI-RS grid. The default value -1 takes the "// &
                          "radius from the extent of the most diffuse AO Gaussian on each atom, "// &
                          "which can be too conservative. A smaller radius can save computation "// &
                          "time, but a too small radius leads to unconverged results.", &
                          usage="CUTOFF_RADIUS_RL_AO 8.0", &
                          default_r_val=cp_unit_to_cp2k(value=-1.0_dp, unit_str="angstrom"), &
                          type_of_var=real_t, unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS_RL_W", &
                          description="Real-space cutoff radius (in Angstrom) for W(r,r') in the "// &
                          "RI-RS GW self-energy Sigma(r,r') = G(r,r') W(r,r'). The Green's function "// &
                          "G(r,r') decays to zero when r and r' are far apart, so W(r,r') is only "// &
                          "needed for |r-r'| below this radius, which saves computation time. "// &
                          "This helps only for very large systems and requires knowing the distance "// &
                          "beyond which G(r,r') is zero. Expert keyword: keep the default -1 "// &
                          "(no cutoff) without such knowledge of the system. GW only: RI_RS_RPA requires -1.", &
                          usage="CUTOFF_RADIUS_RL_W 20.0", &
                          default_r_val=cp_unit_to_cp2k(value=-1.0_dp, unit_str="angstrom"), &
                          type_of_var=real_t, unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_ri_rs_grid_optimization_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_ri_rs_section

! **************************************************************************************************
!> \brief Creates the RI-RS grid optimization section.
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_rs_grid_optimization_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create( &
         section, __LOCATION__, name="GRID_OPTIMIZATION", repeats=.FALSE., &
         description="Optimizes RI-RS grid using small local atomic environments. For every atom "// &
         "in the structure, a cluster is formed from all atoms whose distance from the centre atom "// &
         "is smaller than CUTOFF_ATOMIC_CLUSTER. The grid points are optimized such that the "// &
         "corresponding RI-RS three-centre integrals reproduce the analytic integrals: "// &
         "For a cluster $C_I$, the normalized squared RI-RS error is "// &
         "$E_I = \frac{1}{N_I} \sum_{P\mu\nu \in C_I} \left|(\mu\nu|P) - "// &
         "\sum_{\mathbf{r}_{A\ell} \in C_I} \phi_\mu(\mathbf{r}_{A\ell}) "// &
         "\phi_\nu(\mathbf{r}_{A\ell}) Z_{A\ell P}\right|^2$. The total objective is "// &
         "the average over all local clusters, $E_\mathrm{loc} = \frac{1}{N_\mathrm{at}} "// &
         "\sum_{I=1}^{N_\mathrm{at}} E_I$. Starting points are generated on radial-Lebedev grids "// &
         "and selected by diagonal Cholesky decomposition in independent atomic clusters. "// &
         "If the selection reaches numerical rank before the requested grid size, unused points "// &
         "in the central atom's Voronoi cell are added by geometric maximin distance. "// &
         "Before Voronoi filtering, the candidate grid contains at least "// &
         "max(30, 2 RS_AO_RATIO) times the number of AO functions. RS_AO_RATIO sets the "// &
         "number of grid points per atom. "// &
         "The optimization adapts the grid to the chemical environment and composition of the "// &
         "structure and can permit substantially smaller grids. Compact grids often work well for "// &
         "small basis sets, but RS_AO_RATIO "// &
         "should be tuned and converged for the chosen basis sets and chemical system.")

      CALL keyword_create( &
         keyword, __LOCATION__, name="RS_AO_RATIO", &
         description="Number of real-space grid points $\mathbf{r}_\ell$ per AO basis function. "// &
         "For every atom A, the number of real-space grid points is $N_\mathbf{r}^A=\alpha\, "// &
         "N_\mathrm{AO}^A$, where $\alpha$ is RS_AO_RATIO and $N_\mathrm{AO}^A$ is the number of "// &
         "AO basis functions of atom A. A negative value is the default and sets $\alpha=7$, except "// &
         "for periodic RPA and GW with k-points and the augmented MOLOPT basis sets, where "// &
         "$\alpha$ ranges from 7 to 14: the smaller the basis set, the larger $\alpha$.", &
         usage="RS_AO_RATIO 7", default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUTOFF_ATOMIC_CLUSTER", &
         description="Radius of the atom-centred cluster used to optimize the grid points. "// &
         "All atoms within this distance of the centre atom are included. "// &
         "A negative value selects 1.4 times the distance to the nearest neighbouring atom.", &
         usage="CUTOFF_ATOMIC_CLUSTER -1", type_of_var=real_t, &
         default_r_val=cp_unit_to_cp2k(value=-1.0_dp, unit_str="angstrom"), &
         unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MAX_ITER", &
         description="Maximum number of objective-and-gradient evaluations in the single "// &
         "cooperative Cartesian grid optimization.", &
         usage="MAX_ITER 20", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ri_rs_grid_optimization_section

! **************************************************************************************************
!> \brief Creates the RI-RS RPA section: RPA correlation energy and exact exchange of the SCF
!>        orbitals with periodic RI-RS and the full k-point mesh of DFT%KPOINTS
!> \param section RI_RS_RPA input section
! **************************************************************************************************
   SUBROUTINE create_ri_rs_rpa_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_RS_RPA", &
                          description="RPA total energy of periodic systems with k-point "// &
                          "sampling, computed with the real-space resolution of the identity "// &
                          "(RI-RS): $E^\text{RPA} = E^\text{DFT} - E^\text{xc} + E^\text{EXX} "// &
                          "+ E_\text{c}^\text{RPA}$. The exact exchange energy $E^\text{EXX}$ "// &
                          "and the RPA correlation energy $E_\text{c}^\text{RPA}$ are evaluated "// &
                          "with the SCF orbitals on the k-point mesh of DFT%KPOINTS.", &
                          citations=[Duchemin2019, Duchemin2021], &
                          n_keywords=7, n_subsections=2, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Activates the RI-RS RPA calculation.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL add_minimax_keywords(section)
      CALL add_eps_filter_keyword(section)
      CALL add_ri_keywords(section)

      NULLIFY (subsection)
      CALL create_ri_rs_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_auto_ri_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_ri_rs_rpa_section

END MODULE input_cp2k_ri_rs
