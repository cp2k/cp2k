!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief input sections for the CP2K-native SKALA and GauXC exchange-correlation backends
! **************************************************************************************************
MODULE input_cp2k_skala
   USE input_constants,                 ONLY: do_gpu_auto,&
                                              do_gpu_device,&
                                              do_gpu_host
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
   USE xc_input_constants,              ONLY: skala_gapw_atom_composite_grid,&
                                              skala_gapw_common_grid,&
                                              skala_gapw_cp2k_default,&
                                              skala_gapw_direct_valence,&
                                              skala_gapw_paw_one_center,&
                                              skala_gapw_paw_one_center_split
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_skala'

   PUBLIC :: create_gauxc_section, create_skala_native_section

CONTAINS

! **************************************************************************************************
!> \brief creates the structure of the GAUXC section
!> \param section the section that will be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_gauxc_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="GAUXC", &
                          description="Use exchange-correlation functionals provided by GauXC. "// &
                          "The current CP2K bridge uses GauXC's molecular quadrature. Periodic inputs "// &
                          "are limited to isolated validation cells and are not compact periodic "// &
                          "material calculations.", &
                          n_keywords=25, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL", &
                          description="Functional name used to create the GauXC integrator. "// &
                          "If MODEL is set, this optional base functional defaults to PBE and "// &
                          "does not need to be specified for MODEL SKALA inputs.", &
                          default_c_val="PBE")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL", &
                          description="GauXC Skala model name or path. Use NONE for conventional GauXC. "// &
                          "MODEL SKALA selects the installed Skala-1.1 model when available; "// &
                          "other Skala models are supplied as .fun files or model names. "// &
                          "NLCC pseudopotentials are rejected by the direct molecular GauXC "// &
                          "Skala path.", &
                          usage="MODEL SKALA", &
                          default_c_val="NONE")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
         description="Select the density representation used by an active GauXC model for "// &
         "METHOD GAPW/GAPW_XC kinds with GTH or ECP pseudopotentials. "// &
         "DIRECT_VALENCE treats those kinds as GPW_TYPE and evaluates the model on "// &
         "the direct valence density; this is the default for non-local models such "// &
         "as Skala. PAW_ONE_CENTER combines the smooth and all hard-minus-soft "// &
         "one-center primitive fields before one model evaluation, retaining all 15 pairwise "// &
         "spin-gradient cross terms (14 additional to the smooth-density baseline) and the "// &
         "non-local couplings. PAW_ONE_CENTER_SPLIT is a diagnostic that evaluates smooth, hard, "// &
         "and soft fields separately. CP2K_DEFAULT retains the representation selected "// &
         "by KIND%GPW_TYPE, basis hardness, and DFT%QS%FORCE_PAW. All-electron kinds "// &
         "and conventional GauXC MODEL NONE calculations are unaffected. In a mixed "// &
         "all-electron/pseudopotential calculation, the selector applies per kind: "// &
         "all-electron kinds retain their one-center reconstruction while pseudopotential "// &
         "kinds follow the selected representation. No mixed-system-specific input is needed.", &
         usage="PSEUDOPOTENTIAL_GAPW_REPRESENTATION DIRECT_VALENCE", &
         enum_c_vals=s2a("DIRECT_VALENCE", "PAW_ONE_CENTER", "CP2K_DEFAULT", &
                         "PAW_ONE_CENTER_SPLIT"), &
         enum_i_vals=[skala_gapw_direct_valence, skala_gapw_paw_one_center, &
                      skala_gapw_cp2k_default, skala_gapw_paw_one_center_split], &
         enum_desc=s2a("Use the direct pseudopotential valence density", &
                       "Reconstruct the GAPW primitive fields before model evaluation", &
                       "Use the representation selected by the ordinary CP2K GAPW rules", &
                       "Diagnose the separate smooth plus hard-minus-soft energy expression"), &
         default_i_val=skala_gapw_direct_valence)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID", &
                          description="Atomic grid size for GauXC. Skala defaults to SUPERFINE "// &
                          "unless this keyword is set explicitly.", &
                          usage="GRID fine", &
                          enum_c_vals=s2a("FINE", "ULTRAFINE", "SUPERFINE", "GM3", "GM5"), &
                          enum_i_vals=[1, 2, 3, 4, 5], &
                          enum_desc=s2a("Fine grid (default)", &
                                        "Ultra fine grid", &
                                        "Super fine grid", &
                                        "GM3 grid", &
                                        "GM5 grid"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RADIAL_QUADRATURE", &
                          description="Radial quadrature scheme for GauXC", &
                          usage="RADIAL_QUADRATURE muraknowles", &
                          enum_c_vals=s2a("MURAKNOWLES", "BECKE", "TREUTLERAHLRICHS", "MURRAYHANDYLAMING"), &
                          enum_i_vals=[1, 2, 3, 4], &
                          enum_desc=s2a("Mura-Knowles radial quadrature (default)", &
                                        "Becke radial quadrature", &
                                        "Treutler-Ahlrichs radial quadrature", &
                                        "Murray-Handy-Laming radial quadrature"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRUNING_SCHEME", &
                          description="Pruning scheme for GauXC. Skala defaults to UNPRUNED "// &
                          "unless this keyword is set explicitly.", &
                          usage="PRUNING_SCHEME robust", &
                          enum_c_vals=s2a("ROBUST", "TREUTLER", "UNPRUNED"), &
                          enum_i_vals=[1, 2, 3], &
                          enum_desc=s2a("Robust pruning (default)", &
                                        "Treutler pruning", &
                                        "Unpruned grid"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BATCH_SIZE", &
                          description="Number of quadrature points processed per GauXC "// &
                          "integrator batch. Larger batches reduce launch overhead but increase "// &
                          "device memory use.", &
                          usage="BATCH_SIZE 512", &
                          default_i_val=512)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DEVICE_RUNTIME_FILL_FRACTION", &
                          description="Fraction of currently available GPU memory preallocated by "// &
                          "a GauXC device runtime. Used whenever either LB_EXECUTION_SPACE or "// &
                          "INT_EXECUTION_SPACE is DEVICE. Lower values leave more memory for "// &
                          "Skala Torch inference; higher values can reduce allocator churn.", &
                          usage="DEVICE_RUNTIME_FILL_FRACTION 0.1", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL_ATOM_CHUNK_SIZE", &
                          description="Atom chunk size for GauXC Skala Torch inference. "// &
                          "A positive value evaluates the Torch model atom block by atom block; "// &
                          "zero disables chunking; the default leaves GauXC's model-specific policy "// &
                          "or the GAUXC_ONEDFT_ATOM_CHUNK_SIZE environment variable in control.", &
                          usage="MODEL_ATOM_CHUNK_SIZE 3", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC_REFERENCE", &
                          description="Explicitly opt into the limited isolated-cell GauXC reference "// &
                          "path for periodic CP2K inputs. This path is restricted to Gamma-only, "// &
                          "single-image METHOD GPW calculations with GTH pseudopotentials and uses "// &
                          "GauXC's molecular quadrature; compact periodic materials require a "// &
                          "dedicated periodic GauXC interface.", &
                          usage="PERIODIC_REFERENCE T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLECULAR_VIRIAL", &
                          description="Print a molecular XC virial diagnostic computed from GauXC "// &
                          "nuclear gradients. This is not an analytical periodic stress tensor.", &
                          usage="MOLECULAR_VIRIAL T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLECULAR_VIRIAL_DEBUG", &
                          description="Check the GauXC molecular XC virial diagnostic against a "// &
                          "finite-difference affine scaling of the molecular coordinates.", &
                          usage="MOLECULAR_VIRIAL_DEBUG T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLECULAR_VIRIAL_DEBUG_DX", &
                          description="Dimensionless finite-difference step for MOLECULAR_VIRIAL_DEBUG.", &
                          usage="MOLECULAR_VIRIAL_DEBUG_DX 1.0E-4", &
                          default_r_val=1.0E-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_GPU", &
                          description="Select GPU execution for GauXC. AUTO (default) uses "// &
                          "the GPU when CP2K and GauXC are built with CUDA support. HOST "// &
                          "forces host-only execution. GPU requires a usable GPU and "// &
                          "aborts otherwise. Unless set explicitly, both LB_EXECUTION_SPACE "// &
                          "and INT_EXECUTION_SPACE follow this selection.", &
                          usage="USE_GPU AUTO", &
                          enum_c_vals=s2a("AUTO", "HOST", "GPU"), &
                          enum_i_vals=[do_gpu_auto, do_gpu_host, do_gpu_device], &
                          enum_desc=s2a("Follow the build configuration", &
                                        "Force host execution", &
                                        "Force GPU execution"), &
                          default_i_val=do_gpu_auto)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LB_EXECUTION_SPACE", &
                          description="Execution space for the GauXC load balancer. "// &
                          "Defaults to DEVICE when USE_GPU enables the GPU and HOST "// &
                          "otherwise. An explicit value overrides that default.", &
                          usage="LB_EXECUTION_SPACE host", &
                          enum_c_vals=s2a("HOST", "DEVICE"), &
                          enum_i_vals=[1, 2], &
                          enum_desc=s2a("Run on host", &
                                        "Run on device"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INT_EXECUTION_SPACE", &
                          description="Execution space for the GauXC quadrature integrator. "// &
                          "Defaults to DEVICE when USE_GPU enables the GPU and HOST "// &
                          "otherwise. An explicit value overrides that default.", &
                          usage="INT_EXECUTION_SPACE host", &
                          enum_c_vals=s2a("HOST", "DEVICE"), &
                          enum_i_vals=[1, 2], &
                          enum_desc=s2a("Run on host", &
                                        "Run on device"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LWD_KERNEL", &
                          description="Local work driver kernel used by the GauXC integrator. "// &
                          "AUTO uses SCHEME1-CUTLASS for device SKALA when GauXC provides it "// &
                          "and otherwise uses GauXC's default for the selected execution space. "// &
                          "DEFAULT always uses GauXC's default. "// &
                          "Device builds with CUTLASS support can use SCHEME1-CUTLASS to run "// &
                          "the grouped local potential update through CUTLASS.", &
                          usage="LWD_KERNEL AUTO", &
                          default_c_val="AUTO")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SKALA_RUNTIME", &
                          description="MPI communicator policy for GauXC Skala. AUTO uses "// &
                          "the force-evaluation communicator for closed-shell calculations and a "// &
                          "rank-local replicated runtime for open-shell calculations. MPI forces the "// &
                          "force-evaluation communicator. SELF forces the previous replicated mode.", &
                          usage="SKALA_RUNTIME AUTO", &
                          enum_c_vals=s2a("AUTO", "MPI", "SELF"), &
                          enum_i_vals=[1, 2, 3], &
                          enum_desc=s2a("Closed-shell MPI, open-shell replicated (default)", &
                                        "Use the force-evaluation MPI communicator", &
                                        "Use rank-local replicated GauXC runtimes"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL_GRADIENT_RUNTIME", &
                          description="MPI communicator policy for GauXC Skala nuclear gradients. "// &
                          "AUTO uses distributed analytical gradients when the GauXC PR #222 fix is "// &
                          "verified and energy/VXC uses MPI; otherwise it retains rank-local gradients. "// &
                          "MPI requires verified support and an MPI energy/VXC runtime. SELF forces "// &
                          "rank-local gradients. AUTO never overrides SKALA_RUNTIME's open-shell safeguard.", &
                          usage="MODEL_GRADIENT_RUNTIME AUTO", &
                          enum_c_vals=s2a("AUTO", "MPI", "SELF"), &
                          enum_i_vals=[1, 2, 3], &
                          enum_desc=s2a("Use verified MPI gradients or safe replicated fallback (default)", &
                                        "Use the force-evaluation MPI communicator when supported", &
                                        "Use rank-local replicated gradient runtimes"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FORCE_NUMERICAL_GRADIENT", &
                          description="Force central finite differences for the GauXC molecular XC "// &
                          "nuclear gradient, using a displacement of 1.0E-4 bohr at fixed AO density. "// &
                          "Diagnostic option; costs six XC energy evaluations per atom. Does not "// &
                          "affect SKALA_NATIVE or the CP2K one-center contribution.", &
                          usage="FORCE_NUMERICAL_GRADIENT T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUTPUT_PATH", &
                          description="Optional path to an existing directory for GauXC HDF5 debug output. "// &
                          "If set, molecule and basis set data are written as separate HDF5 files.", &
                          usage="OUTPUT_PATH /path/to/output", &
                          default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_gauxc_section

! **************************************************************************************************
!> \brief creates the structure of the CP2K-native SKALA section
!> \param section the section that will be created
! **************************************************************************************************
   SUBROUTINE create_skala_native_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SKALA_NATIVE", &
                          description="Evaluate the SKALA TorchScript model on CP2K's native "// &
                          "real-space grid. This is the CP2K-native SKALA/FTorch implementation "// &
                          "and is mutually exclusive with the GAUXC section. It supports METHOD "// &
                          "GPW, GAPW, and GAPW_XC, including analytical gradients and stress.", &
                          n_keywords=25, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the CP2K-native SKALA implementation", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL", &
                          description="SKALA TorchScript model for the CP2K-native grid path. "// &
                          "MODEL SKALA resolves the installed SKALA model; an explicit path "// &
                          "selects a .fun file. With GPU execution the model is resolved from "// &
                          "SKALA_CUDA_MODEL, otherwise from SKALA_MODEL.", &
                          usage="MODEL SKALA", &
                          default_c_val="SKALA")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
         description="Select the density representation used by the CP2K-native SKALA model for "// &
         "METHOD GAPW/GAPW_XC kinds with GTH or ECP pseudopotentials. "// &
         "DIRECT_VALENCE treats those kinds as GPW_TYPE and evaluates the model on "// &
         "the direct valence density; this is the default for non-local models such "// &
         "as Skala. PAW_ONE_CENTER combines the smooth and all hard-minus-soft "// &
         "one-center primitive fields before one model evaluation, retaining all 15 pairwise "// &
         "spin-gradient cross terms (14 additional to the smooth-density baseline) and the "// &
         "non-local couplings. PAW_ONE_CENTER_SPLIT is a diagnostic that evaluates smooth, hard, "// &
         "and soft fields separately. CP2K_DEFAULT retains the representation selected "// &
         "by KIND%GPW_TYPE, basis hardness, and DFT%QS%FORCE_PAW. All-electron kinds "// &
         "are unaffected. In a mixed "// &
         "all-electron/pseudopotential calculation, the selector applies per kind: "// &
         "all-electron kinds retain their one-center reconstruction while pseudopotential "// &
         "kinds follow the selected representation. No mixed-system-specific input is needed.", &
         usage="PSEUDOPOTENTIAL_GAPW_REPRESENTATION DIRECT_VALENCE", &
         enum_c_vals=s2a("DIRECT_VALENCE", "PAW_ONE_CENTER", "CP2K_DEFAULT", &
                         "PAW_ONE_CENTER_SPLIT"), &
         enum_i_vals=[skala_gapw_direct_valence, skala_gapw_paw_one_center, &
                      skala_gapw_cp2k_default, skala_gapw_paw_one_center_split], &
         enum_desc=s2a("Use the direct pseudopotential valence density", &
                       "Reconstruct the GAPW primitive fields before model evaluation", &
                       "Use the representation selected by the ordinary CP2K GAPW rules", &
                       "Diagnose the separate smooth plus hard-minus-soft energy expression"), &
         default_i_val=skala_gapw_direct_valence)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL_FP64", &
                          description="Promote float32 numerical state and explicit floating-point casts "// &
                          "in the native SKALA TorchScript model to double precision. Trained parameter "// &
                          "values are retained exactly. This improves SCF and geometry-derivative "// &
                          "consistency for mixed-precision exports. Disable to retain the export's "// &
                          "original numerical precision.", &
                          usage="MODEL_FP64 T", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_GPU", &
                          description="Select GPU execution for the CP2K-native SKALA "// &
                          "TorchScript model. AUTO (default) uses the GPU when CP2K is "// &
                          "built with a CUDA-enabled LibTorch. HOST forces CPU execution. "// &
                          "GPU requires a usable GPU and aborts otherwise. GPU execution "// &
                          "resolves MODEL SKALA from SKALA_CUDA_MODEL; HOST uses SKALA_MODEL.", &
                          usage="USE_GPU AUTO", &
                          enum_c_vals=s2a("AUTO", "HOST", "GPU"), &
                          enum_i_vals=[do_gpu_auto, do_gpu_host, do_gpu_device], &
                          enum_desc=s2a("Follow the build configuration", &
                                        "Force CPU execution", &
                                        "Force GPU execution"), &
                          default_i_val=do_gpu_auto)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_PARTITION", &
                          description="Atom partition used to build the experimental CP2K-native "// &
                          "SKALA COMMON_GRID feature block. HARD assigns each native grid point to "// &
                          "the nearest periodic atom and is intended for legacy energy/VXC-only runs. "// &
                          "SMOOTH uses a Becke-like fuzzy-cell partition on the native grid and creates "// &
                          "weighted atom rows. Internal atom-grid weights vanish through a smooth taper "// &
                          "when sparse rows reach the partition cutoff. Analytical atom forces and stress "// &
                          "include both the partition and taper derivatives in the SKALA feature response.", &
                          usage="ATOM_PARTITION SMOOTH", &
                          enum_c_vals=s2a("HARD", "SMOOTH"), &
                          enum_i_vals=[1, 2], &
                          enum_desc=s2a("Assign each grid point to the nearest atom", &
                                        "Use a smooth Becke-like native-grid atom partition"), &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="LAYOUT", &
         description="Select the numerical layout for native-grid Skala irrespective of the "// &
         "electronic core representation. ATOM_COMPOSITE, the default, evaluates complete "// &
         "atom-centered radial/Lebedev blocks. GPW and pseudopotential direct-valence fields "// &
         "are interpolated from CP2K's plane-wave grid; GAPW primitive fields are combined "// &
         "before nonlinear feature construction. The matching interpolation, partition, and "// &
         "one-center adjoints return VXC, forces, and stress to CP2K. Mixed all-electron and "// &
         "pseudopotential kinds use the same syntax and are reconstructed per kind. COMMON_GRID retains "// &
         "evaluation on CP2K's regular real-space grid as a cutoff-sensitive reference and "// &
         "diagnostic layout; mixed core representations require ATOM_COMPOSITE.", &
         usage="LAYOUT ATOM_COMPOSITE", &
         enum_c_vals=s2a("COMMON_GRID", "ATOM_COMPOSITE"), &
         enum_i_vals=[skala_gapw_common_grid, skala_gapw_atom_composite_grid], &
         enum_desc=s2a("Evaluate on the common regular real-space grid", &
                       "Use image-partitioned atom-centered grids"), &
         default_i_val=skala_gapw_atom_composite_grid)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_DENSITY_PARTITION", &
                          description="Hard/soft one-center primitive-field contribution used by "// &
                          "METHOD GAPW and METHOD GAPW_XC. With PAW_ONE_CENTER, HARD_MINUS_SOFT "// &
                          "combines smooth, hard, and soft density, gradient, and kinetic-energy-density "// &
                          "fields before a single nonlinear SKALA evaluation. It is the default. For "// &
                          "pseudopotential PAW/one-center GAPW this is distinct from KIND%GPW_TYPE, "// &
                          "which uses the regular-grid valence-density approximation. "// &
                          "HARD_ONLY, SOFT_ONLY, and NONE are diagnostic variants.", &
                          usage="GAPW_DENSITY_PARTITION HARD_MINUS_SOFT", &
                          enum_c_vals=s2a("HARD_MINUS_SOFT", "HARD_ONLY", "SOFT_ONLY", "NONE"), &
                          enum_i_vals=[1, 2, 3, 4], &
                          enum_desc=s2a("Use the hard-minus-soft GAPW one-center correction", &
                                        "Use only the hard one-center density", &
                                        "Use only the soft one-center density with the GAPW sign", &
                                        "Disable the GAPW one-center SKALA correction"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_COMPOSITE_REFERENCE", &
                          description="Evaluate native-grid SKALA for METHOD GAPW or GAPW_XC on "// &
                          "rho and tau reconstructed as smooth plus hard-minus-soft before nonlinear "// &
                          "feature construction on one common regular grid. No separate one-center "// &
                          "SKALA energy is added. This expensive, cutoff-sensitive diagnostic is "// &
                          "intended to validate the nonlinear hard/soft cross terms.", &
                          usage="GAPW_COMPOSITE_REFERENCE T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_COMPOSITE_DIRECT_AO", &
                          description="Use direct full-ORB rho and tau collocation for the native-grid "// &
                          "GAPW composite diagnostic. This is the variational all-electron "// &
                          "counterpart of the reconstructed smooth plus hard-minus-soft fields and "// &
                          "provides a reference for their common-grid adjoint.", &
                          usage="GAPW_COMPOSITE_DIRECT_AO T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_ATOM_COMPOSITE_REFERENCE", &
                          description="Replace the native smooth plus separate GAPW one-center SKALA "// &
                          "terms by one atom-centered evaluation of smooth plus hard-minus-soft "// &
                          "primitive fields before nonlinear feature construction. The smooth PW "// &
                          "and one-center VXC adjoints are returned consistently. MPI ranks assemble "// &
                          "the atom rows collectively; one communicator root evaluates the model and "// &
                          "broadcasts the energy and feature adjoints.", &
                          usage="GAPW_ATOM_COMPOSITE_REFERENCE T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_ATOM_COMPOSITE_DIAGNOSTIC", &
                          description="Evaluate an energy-only GAPW diagnostic on the combined "// &
                          "smooth plus one-center hard-minus-soft density before nonlinear SKALA "// &
                          "feature construction. The smooth native-grid density is interpolated to "// &
                          "the atom-centered GAPW quadrature and all atomic grids are evaluated in "// &
                          "one model call. This single-rank diagnostic does not replace the active "// &
                          "XC energy or potential.", &
                          usage="GAPW_ATOM_COMPOSITE_DIAGNOSTIC T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="GAPW_ATOM_COMPOSITE_COMPONENTS", &
                          description="Select which one-center hard-minus-soft fields enter the "// &
                          "atom-composite SKALA diagnostic. FULL combines rho, gradient, and tau; "// &
                          "RHO_GRAD and TAU isolate the corresponding reconstruction; SMOOTH "// &
                          "uses only interpolated smooth fields.", &
                          usage="GAPW_ATOM_COMPOSITE_COMPONENTS FULL", &
                          enum_c_vals=s2a("FULL", "RHO_GRAD", "TAU", "SMOOTH"), &
                          enum_i_vals=[1, 2, 3, 4], &
                          enum_desc=s2a("Combine rho, gradient, and tau corrections", &
                                        "Combine only rho and gradient corrections", &
                                        "Combine only the tau correction", &
                                        "Use only interpolated smooth fields"), &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="GAPW_ATOM_COMPOSITE_FD_STEP", &
                          description="Central finite-difference step used by the atom-composite "// &
                          "feature-VXC diagnostic.", &
                          usage="GAPW_ATOM_COMPOSITE_FD_STEP 3.0E-3", &
                          default_r_val=3.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_CHUNKS", &
                          description="Split the CP2K-native SKALA TorchScript "// &
                          "evaluation into contiguous atom blocks. MPI runs evaluate one atom block per "// &
                          "rank, while single-rank runs can still use atom subchunks to reduce peak CUDA "// &
                          "memory. Dynamic feature rows and VXC gradients are routed between their "// &
                          "local grid and atom-block owner ranks. This is enabled by default; set it to "// &
                          "false to retain the replicated global feature block. ATOM_COMPOSITE also "// &
                          "subchunks energy, VXC, force, and stress evaluations along independent atomic "// &
                          "blocks; COMMON_GRID analytical force and stress runs retain the full "// &
                          "differentiable tensor graph.", &
                          usage="ATOM_CHUNKS T", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_CHUNK_MAX_ROWS", &
                          description="Maximum number of atom-grid rows evaluated by one TorchScript "// &
                          "call in the experimental CP2K-native SKALA GPW atom-chunk path. A positive "// &
                          "value splits each rank-local atom chunk into contiguous atom subchunks "// &
                          "whose padded Torch row count stays below the limit. Atomic blocks are not "// &
                          "divided, so one block can exceed the requested limit. Zero disables "// &
                          "subchunking; -1 selects an automatic padded-row cap.", &
                          usage="ATOM_CHUNK_MAX_ROWS 250000", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIAGNOSTICS", &
                          description="Print diagnostic integrals of the CP2K-native SKALA GPW feature "// &
                          "block, including electron count, spin moment, and summed grid weights.", &
                          usage="DIAGNOSTICS T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_skala_native_section

END MODULE input_cp2k_skala
