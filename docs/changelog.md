# Changelog

## 2026.2 (July 15, 2026)

### New Features

- DFT+U with k-points for Mulliken methods ([#4855](https://github.com/cp2k/cp2k/pull/4855))
- Energy Correction Harris functional with k-points
  ([#5031](https://github.com/cp2k/cp2k/pull/5031))
- Lowdin population analysis for k-points ([#5045](https://github.com/cp2k/cp2k/pull/5045))
- Wavefunction extrapolation for k-point calculations
  ([#4884](https://github.com/cp2k/cp2k/pull/4884), [#4943](https://github.com/cp2k/cp2k/pull/4943),
  [#4949](https://github.com/cp2k/cp2k/pull/4949), [#5219](https://github.com/cp2k/cp2k/pull/5219))
- GExt wavefunction extrapolation ([#5043](https://github.com/cp2k/cp2k/pull/5043),
  [#5229](https://github.com/cp2k/cp2k/pull/5229))
- Adaptively Compressed Exchange (ACE) option to the HFX/ADMM ground-state path
  ([#5238](https://github.com/cp2k/cp2k/pull/5238))
- Alternative smearing methods ([#4958](https://github.com/cp2k/cp2k/pull/4958))
- Brownian chain molecular dynamics (BCMD) propagator for path integrals
  ([#5247](https://github.com/cp2k/cp2k/pull/5247))
- DOS/PDOS with broadened output and k-point projections
  ([#5287](https://github.com/cp2k/cp2k/pull/5287),
  [#5299](https://github.com/cp2k/cp2k/pull/5299),)
- K-point symmetry reduction ([#5123](https://github.com/cp2k/cp2k/pull/5123),
  [#5152](https://github.com/cp2k/cp2k/pull/5152), [#5165](https://github.com/cp2k/cp2k/pull/5165),
  [#5173](https://github.com/cp2k/cp2k/pull/5173))
- Input keyword for cutoff radius of HFX shortrange potential
  ([#4945](https://github.com/cp2k/cp2k/pull/4945))
- Input keyword for L-BFGS optimizer print control ([#5274](https://github.com/cp2k/cp2k/pull/5274))
- New output of NEB including relative energy plot and final structures
  ([#5382](https://github.com/cp2k/cp2k/pull/5382))
- New output of final structure from optimization in CIF and EXTXYZ formats
  ([#5118](https://github.com/cp2k/cp2k/pull/5118), [#5140](https://github.com/cp2k/cp2k/pull/5140))
- New output of Hessian before mass weighting from revised vibrational analysis printout
  ([#5395](https://github.com/cp2k/cp2k/pull/5395))
- Dimer initialization by reading molden files ([#5312](https://github.com/cp2k/cp2k/pull/5312))
- Parse cell information from EXTXYZ for both initial geometry and each frame of reftraj run
  ([#4960](https://github.com/cp2k/cp2k/pull/4960), [#5401](https://github.com/cp2k/cp2k/pull/5401))
- Wrapping input atomic coordinates for each frame before calculation in a reftraj run
  ([#5479](https://github.com/cp2k/cp2k/pull/5479))
- Per-thermal-region function for rescaling temperatures in MD (independent of thermostats)
  ([#5002](https://github.com/cp2k/cp2k/pull/5002))
- Cell optimization with fixed volume ([#5086](https://github.com/cp2k/cp2k/pull/5086))

### New Libraries

- Add openPMD output, including scientific metadata support
  ([#4058](https://github.com/cp2k/cp2k/pull/4058), [#4931](https://github.com/cp2k/cp2k/pull/4931),
  [#4942](https://github.com/cp2k/cp2k/pull/4942), [#5096](https://github.com/cp2k/cp2k/pull/5096))
- GauXC/Skala models ([#5084](https://github.com/cp2k/cp2k/pull/5084))
- LibFCI active-space solver ([#5167](https://github.com/cp2k/cp2k/pull/5167))
- Reintegrate with LIBXS, LIBXSTREAM, and LIBXSMM ([#5343](https://github.com/cp2k/cp2k/pull/5343))
- Add libGint for Hartree–Fock exchange with CUDA acceleration
  ([#5446](https://github.com/cp2k/cp2k/pull/5446))

### Breaking Changes

- Remove &MOTION/&CELL_OPT/TYPE; cell optimizations now always use DIRECT_CELL_OPT
  ([#5257](https://github.com/cp2k/cp2k/pull/5257))
- Drop support for GCC 8 ([#5290](https://github.com/cp2k/cp2k/pull/5290))
- Refactor DOS/PDOS input section ([#5326](https://github.com/cp2k/cp2k/pull/5326))
- Remove obsolete Cython Python bindings ([#5541](https://github.com/cp2k/cp2k/pull/5541))
- Remove EVAL_ENERGY_FORCES and EVAL_FORCES keywords in favor of EVAL under &MOTION/&MD/&REFTRAJ
  ([#5401](https://github.com/cp2k/cp2k/pull/5401))
- An implementation of the FFTW3 interface may be turned into a hard dependency in a later release.
  Please consider compiling and CP2K with FFTW3, MKL, AOCL or any other library implementing this
  interface if you have not used it until now. ([#5454](https://github.com/cp2k/cp2k/pull/5454))

### Fixes

- Fix issues in the native DFT-D4 implementation. ([#5030](https://github.com/cp2k/cp2k/pull/5030))
- Analytical periodic-subspace stress for 2D systems with ANALYTIC/MT
  ([#5282](https://github.com/cp2k/cp2k/pull/5282))

______________________________________________________________________

## 2026.1 (January 6, 2026)

### New Features

- Add new BASIS_AUG_MOLOPT and BASIS_RI_AUG_MOLOPT basis sets (all electron)
  ([#4354](https://github.com/cp2k/cp2k/pull/4354))
- Active Space: Added a new longrange truncated potential for the ERI
  ([#4357](https://github.com/cp2k/cp2k/pull/4357))
- Activate handling of GTH-SOC PP in Sirius ([#4363](https://github.com/cp2k/cp2k/pull/4363))
- Better availability of occupied orbital eigenvalues from OT
  ([#4364](https://github.com/cp2k/cp2k/pull/4364))
- CNEO-DFT energy and force ([#4403](https://github.com/cp2k/cp2k/pull/4403),
  [#4420](https://github.com/cp2k/cp2k/pull/4420))
- Sternheimer BSE (1. Version with W matrix given) ([#4445](https://github.com/cp2k/cp2k/pull/4445))
- SF-TDDFT implementation ([#4446](https://github.com/cp2k/cp2k/pull/4446))
- PAO-ML: Migrate to NequIP framework ([#4461](https://github.com/cp2k/cp2k/pull/4461))
- GAPW DC-DFT + EC EXTERNAL ([#4502](https://github.com/cp2k/cp2k/pull/4502))
- NEGF: New method to extract the matrix Hamiltonians for electrodes
  ([#4520](https://github.com/cp2k/cp2k/pull/4520), [#4636](https://github.com/cp2k/cp2k/pull/4636))
- precommit: Add Fortitude linter ([#4523](https://github.com/cp2k/cp2k/pull/4523))
- RT-TDDFT+RT-BSE : Fourier Transform Outputs ([#4527](https://github.com/cp2k/cp2k/pull/4527))
- RIXS: Open-Shell ([#4533](https://github.com/cp2k/cp2k/pull/4533))
- Add interface to the MiMiC framework for multiscale simulations
  ([#4546](https://github.com/cp2k/cp2k/pull/4546))
- Restricted Space Excitations for TDDFPT ([#4560](https://github.com/cp2k/cp2k/pull/4560),
  [#4588](https://github.com/cp2k/cp2k/pull/4588))
- Vibronic spectroscopy post-processing tools ([#4581](https://github.com/cp2k/cp2k/pull/4581))
- Add extended XYZ format for trajectories ([#4601](https://github.com/cp2k/cp2k/pull/4601))
- Add MAX_FILE_SIZE_MB keyword to MO_CUBES section ([#4603](https://github.com/cp2k/cp2k/pull/4603))
- Implementation of moments for k-point calculations
  ([#4621](https://github.com/cp2k/cp2k/pull/4621))

### Breaking Changes

- Remove Makefile ([#4618](https://github.com/cp2k/cp2k/pull/4618))
- Remove QUIP ([#4616](https://github.com/cp2k/cp2k/pull/4616))
- Some remaining issues with Intel compilers and CMake
  ([#4550](https://github.com/cp2k/cp2k/pull/4550))

### Fixes

- Fix missing (zero) orbital eigenvalues in the MOLDEN output
  ([#4364](https://github.com/cp2k/cp2k/pull/4364))
- Fix gradient for peridic systems with tblite([#4397](https://github.com/cp2k/cp2k/pull/4397))
- Fix stress tensor in tblite ([#4406](https://github.com/cp2k/cp2k/pull/4406))
- Fix force constant in vibrational analysis ([#4427](https://github.com/cp2k/cp2k/pull/4427))
- Fix molecule indexing when using MULTIPLE_UNIT_CELL
  ([#4437](https://github.com/cp2k/cp2k/pull/4437))
- Fix HFX initialization performance issue ([#4455](https://github.com/cp2k/cp2k/pull/4455))
- Fix 3PNT linesearch for CG in GEO_OPT ([#4665](https://github.com/cp2k/cp2k/pull/4665))

______________________________________________________________________

## 2025.2 (July 23, 2025)

### New Features

- TDDFPT: Add exciton descriptors ([#3847](https://github.com/cp2k/cp2k/pull/3847))
- Efield implementation for xTB ([#3883](https://github.com/cp2k/cp2k/pull/3883))
- Atom code: Add option to write xmgrace wavefunction file
  ([#3962](https://github.com/cp2k/cp2k/pull/3962))
- GFN-xTB ([#4005](https://github.com/cp2k/cp2k/pull/4005),
  [#4190](https://github.com/cp2k/cp2k/pull/4190), [#4209](https://github.com/cp2k/cp2k/pull/4209),
  [#4242](https://github.com/cp2k/cp2k/pull/4242))
- Read/write external energy derivatives from/to trexio files
  ([#4009](https://github.com/cp2k/cp2k/pull/4009), [#4074](https://github.com/cp2k/cp2k/pull/4074))
- Pseudopotentials: Add 5f-in-core for trivalent and tetravalent actinides
  ([#4068](https://github.com/cp2k/cp2k/pull/4068))
- Pseudopotentials: Add ccECP ([#3940](https://github.com/cp2k/cp2k/pull/3940))
- RTBSE : Padé FT Refinement ([#4115](https://github.com/cp2k/cp2k/pull/4115))
- HP-DFT modules and regtests ([#4138](https://github.com/cp2k/cp2k/pull/4138))
- Add option to print space groups ([#4271](https://github.com/cp2k/cp2k/pull/4271))
- Add SIRIUS DFTD3 and DFTD4 support ([#4277](https://github.com/cp2k/cp2k/pull/4277))
- Atomic polarization tensors via numerical differentiation
  ([#4287](https://github.com/cp2k/cp2k/pull/4287))
- RI-HFXk: Various improvements ([#4291](https://github.com/cp2k/cp2k/pull/4291))
- Resonant Inelastic X-ray Spectroscopy (RIXS) Module
  ([#4315](https://github.com/cp2k/cp2k/pull/4315))

### New Libraries

- Improve DLA-Future integration ([#4169](https://github.com/cp2k/cp2k/pull/4169),
  [#4269](https://github.com/cp2k/cp2k/pull/4269))
- Upgrade to DeePMD 3.1.0 and switch to PyTorch backend
  ([#3893](https://github.com/cp2k/cp2k/pull/3893), [#4310](https://github.com/cp2k/cp2k/pull/4310))
- Make GRPP an internal dependency ([#3966](https://github.com/cp2k/cp2k/pull/3966))
- Interface for greenX library ([#4078](https://github.com/cp2k/cp2k/pull/4078))
- Add ace support ([#4182](https://github.com/cp2k/cp2k/pull/4182))

### Breaking Changes

- Remove old TDDFPT code ([#4066](https://github.com/cp2k/cp2k/pull/4066))
- Restore old format for writing forces to .xyz files
  ([#4294](https://github.com/cp2k/cp2k/pull/4294))
- RTBSE: Input structure changed ([#3918](https://github.com/cp2k/cp2k/pull/3918))

### Fixes

- Fix occurrence of ghost states with SVDs ([#3911](https://github.com/cp2k/cp2k/pull/3911))
- Speedup loading of NequIP and Allegro models ([#3912](https://github.com/cp2k/cp2k/pull/3912))
- TDDFPT%MGRID bug fix for GAPW ([#3967](https://github.com/cp2k/cp2k/pull/3967))
- Fix typo in data/BASIS_ccGRB_UZH ([#4017](https://github.com/cp2k/cp2k/pull/4017))

______________________________________________________________________

## 2025.1 (January 1, 2025)

### New Features

- BSE: Optical spectra, BSE@evGW(0), and lots more ([#3628](https://github.com/cp2k/cp2k/pull/3628),
  [#3659](https://github.com/cp2k/cp2k/pull/3659), [#3781](https://github.com/cp2k/cp2k/pull/3781),
  [#3793](https://github.com/cp2k/cp2k/pull/3793), [#3815](https://github.com/cp2k/cp2k/pull/3815))
- Real time Bethe-Salpeter Propagation ([#3691](https://github.com/cp2k/cp2k/pull/3691))
- Harris and EHT methods incl. LS solver ([#3665](https://github.com/cp2k/cp2k/pull/3665),
  [#3780](https://github.com/cp2k/cp2k/pull/3780), [#3790](https://github.com/cp2k/cp2k/pull/3790))
- Print wannier states coefficients in AO basis ([#3683](https://github.com/cp2k/cp2k/pull/3683),
  [#3687](https://github.com/cp2k/cp2k/pull/3687))
- Add Z-matrix formalism and EVERY_N_STEP keyword for linear response
  ([#3689](https://github.com/cp2k/cp2k/pull/3689), [#3692](https://github.com/cp2k/cp2k/pull/3692))
- RPA based SIGMA Functionals ([#3695](https://github.com/cp2k/cp2k/pull/3695))
- Framework to calculate response forces for external energy expressions based on KS orbitals
  ([#3721](https://github.com/cp2k/cp2k/pull/3721))
- PAO: Add prediction from equivariant PyTorch models
  ([#3738](https://github.com/cp2k/cp2k/pull/3738))
- gfn0-xTB and parallel DFT-D4 ([#3765](https://github.com/cp2k/cp2k/pull/3765),
  [#3679](https://github.com/cp2k/cp2k/pull/3679), [#3685](https://github.com/cp2k/cp2k/pull/3685))
- Smeared occupation TDA ([#3829](https://github.com/cp2k/cp2k/pull/3829))
- GW: Add keywords SIZE_LATTICE_SUM and KPOINTS_W ([#3833](https://github.com/cp2k/cp2k/pull/3833))

### New Libraries

- Add [SMEAGOL](https://github.com/StefanoSanvitoGroup/libsmeagol) for electron transport with NEGF
  ([#3716](https://github.com/cp2k/cp2k/pull/3716))
- Add [cuSOLVERMp](https://manual.cp2k.org/trunk/technologies/eigensolvers/cusolvermp.html)
  generalized eigensolver ([#3787](https://github.com/cp2k/cp2k/pull/3787))
- Add [DLA-Future](https://manual.cp2k.org/trunk/technologies/eigensolvers/dlaf.html) generalized
  and complex eigensolvers ([#3799](https://github.com/cp2k/cp2k/pull/3799),
  [#3813](https://github.com/cp2k/cp2k/pull/3813), [#3819](https://github.com/cp2k/cp2k/pull/3819))
- Add [TREXIO](https://github.com/TREX-CoE/trexio) file format writer
  ([#3792](https://github.com/cp2k/cp2k/pull/3792))

### Breaking Changes

- Weighting of RSMD colvar is modified for the case of subsystem = list
  ([#3818](https://github.com/cp2k/cp2k/pull/3818))

### Fixes

- Re-init FFT scratch in case of PW env changes ([#3661](https://github.com/cp2k/cp2k/pull/3661))
- Add asserts after malloc and realloc in DBM and grid code
  ([#3706](https://github.com/cp2k/cp2k/pull/3706))
- Fix segmentation fault when performing CDFT-MD with many constraints
  ([#3711](https://github.com/cp2k/cp2k/pull/3711))
- Fix thread-safety in fft_tools ([#3729](https://github.com/cp2k/cp2k/pull/3729))
- Fix access to unallocated arrays in pme, spme, and ewalds
  ([#3733](https://github.com/cp2k/cp2k/pull/3733))
- Fix handling of empty file names in the input ([#3758](https://github.com/cp2k/cp2k/pull/3758))
- Fix bug in TDDFPT/Davidson restart ([#3821](https://github.com/cp2k/cp2k/pull/3821))
- Fix bug in RTP/EMD restart ([#3637](https://github.com/cp2k/cp2k/pull/3637))
- Add initialization step to i-Pi which some clients expect
  ([#3747](https://github.com/cp2k/cp2k/pull/3747))

______________________________________________________________________

## 2024.3 (September 9, 2024)

This is a minor release to fix an issue with the PW environment that can lead to stalls during MD
(https://github.com/cp2k/cp2k/issues/3661).

Since this only affects MPI jobs, the ssmp binaries from the previous
[2024.2 release](https://github.com/cp2k/cp2k/releases/tag/v2024.2) are still up-to-date.

______________________________________________________________________

## 2024.2 (August 6, 2024)

### New Features

- Many new pages in the methods section of [the manual](https://manual.cp2k.org/trunk/)
- ECP nuclear gradients ([#3210](https://github.com/cp2k/cp2k/pull/3210))
- New basis sets ([#3266](https://github.com/cp2k/cp2k/pull/3266),
  [#3276](https://github.com/cp2k/cp2k/pull/3276) [#3460](https://github.com/cp2k/cp2k/pull/3460),
  [#3561](https://github.com/cp2k/cp2k/pull/3561))
- TDA Kernel methods ([#3273](https://github.com/cp2k/cp2k/pull/3273))
- Bethe Salpeter equation for molecules ([#3308](https://github.com/cp2k/cp2k/pull/3308),
  [#3329](https://github.com/cp2k/cp2k/pull/3329))
- Stress prediction in NequIP & Allegro ([#3428](https://github.com/cp2k/cp2k/pull/3428),
  [#3445](https://github.com/cp2k/cp2k/pull/3445))
- GW for small unit cells with full k-points ([#3448](https://github.com/cp2k/cp2k/pull/3448),
  [#3505](https://github.com/cp2k/cp2k/pull/3505), [#3513](https://github.com/cp2k/cp2k/pull/3513))
- Hedin shift for G0W0 ([#3533](https://github.com/cp2k/cp2k/pull/3533))
- i-PI server functionality ([#3420](https://github.com/cp2k/cp2k/pull/3420))
- Infrastructure for Kpoint symmetries ([#3482](https://github.com/cp2k/cp2k/pull/3482))

### New Libraries

- Machine learning with the [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit)
  ([#3145](https://github.com/cp2k/cp2k/pull/3145))
- Dispersion correction with the [DFTD4 library](https://github.com/dftd4/dftd4)
  ([#3501](https://github.com/cp2k/cp2k/pull/3501))
- GPU support via OpenCL ([#3321](https://github.com/cp2k/cp2k/pull/3321),
  [#3315](https://github.com/cp2k/cp2k/pull/3315), [#3375](https://github.com/cp2k/cp2k/pull/3375))

### Breaking Changes

- Increase ScaLAPACK default block size to 64 ([#3184](https://github.com/cp2k/cp2k/pull/3184))
- Remove `BROYDEN_MIXING_NEW` option ([#3346](https://github.com/cp2k/cp2k/pull/3346))
- Remove `KP_RI_EXTENSION_FACTOR"` keyword ([#3223](https://github.com/cp2k/cp2k/pull/3223))
- Mark support for QUIP and PEXSI as deprecated ([#3600](https://github.com/cp2k/cp2k/pull/3600))

### Fixes

- Fix oscillator strength for TDDFT+SOC within length and velocity representation
  ([#3201](https://github.com/cp2k/cp2k/pull/3201))
- External potential: Do not re-parse potential function at each grid-point
  ([#3204](https://github.com/cp2k/cp2k/pull/3204))
- PWDFT: Print the band gap and total energy ([#3211](https://github.com/cp2k/cp2k/pull/3211))
- Fix bug [#3218](https://github.com/cp2k/cp2k/pull/3218) in screening of hfx derivatives
  ([#3221](https://github.com/cp2k/cp2k/pull/3221))
- Fix bug [#3217](https://github.com/cp2k/cp2k/pull/3217) in printout of eigenvalues and
  eigenvectors ([#3230](https://github.com/cp2k/cp2k/pull/3230))
- Improve element comparison in reftraj ([#3337](https://github.com/cp2k/cp2k/pull/3337))
- Allegro: Fix parallelization with virials ([#3445](https://github.com/cp2k/cp2k/pull/3445))
- FFTW: Fix import/export wisdom file ([#3492](https://github.com/cp2k/cp2k/pull/3492),
  [#3496](https://github.com/cp2k/cp2k/pull/3496))

______________________________________________________________________

## 2024.1 (January 3, 2024)

### New Features

- Docs: Launch Sphinx-based [manual](https://manual.cp2k.org)
  ([#2883](https://github.com/cp2k/cp2k/pull/2883))
- TDDFPT: Pseudopotentials, GAPW, GPW, and forces ([#2895](https://github.com/cp2k/cp2k/pull/2895),
  [#2932](https://github.com/cp2k/cp2k/pull/2932) , [#2940](https://github.com/cp2k/cp2k/pull/2940),
  [#3011](https://github.com/cp2k/cp2k/pull/3011), [#3110](https://github.com/cp2k/cp2k/pull/3110),
  [#3122](https://github.com/cp2k/cp2k/pull/3122))
- RI-HFX for K-points (with gradients and ADMM) ([#2998](https://github.com/cp2k/cp2k/pull/2998))
- ADMM: input short cuts ([#3118](https://github.com/cp2k/cp2k/pull/3118))
- Long-range quantum computing (WF) in short-range DFT embedding
  ([#2924](https://github.com/cp2k/cp2k/pull/2924))
- Active space: Implement ERI calculation using half-transformed integrals
  ([#3082](https://github.com/cp2k/cp2k/pull/3082))
- GW: open-shell periodic GW ([#2920](https://github.com/cp2k/cp2k/pull/2920),
  [#3045](https://github.com/cp2k/cp2k/pull/3045), [#3125](https://github.com/cp2k/cp2k/pull/3125))
- G0W0/SOC: bandstructure, PDOS, local bandgap, periodic
  ([#2994](https://github.com/cp2k/cp2k/pull/2994), [#3130](https://github.com/cp2k/cp2k/pull/3130))
- NNP: Helium-Solute for interaction ([#3043](https://github.com/cp2k/cp2k/pull/3043))
- Time-dependent Field for EMD ([#3081](https://github.com/cp2k/cp2k/pull/3081))

### New Libraries

- Add experimental support for [DLA-Future](https://github.com/eth-cscs/DLA-Future) eigensolver
  ([#3143](https://github.com/cp2k/cp2k/pull/3143))
- Enabling calculations with ECPs [libgrpp](https://github.com/aoleynichenko/libgrpp) library
  ([#3147](https://github.com/cp2k/cp2k/pull/3147))

### Breaking Changes

- Remove SINGLE_PRECISION_MATRICES keyword ([#3096](https://github.com/cp2k/cp2k/pull/3096),
  [#3140](https://github.com/cp2k/cp2k/pull/3140))
- Abort run by default on SCF convergence failure ([#3148](https://github.com/cp2k/cp2k/pull/3148))
- Drop support for NDEBUG ([#3172](https://github.com/cp2k/cp2k/pull/3172))
- Production docker files moved to [new repository](https://github.com/cp2k/cp2k-containers)
  ([#3083](https://github.com/cp2k/cp2k/pull/3083))
- MD: Refactor REFTRAJ / EVAL_ENERGY_FORCES keyword
  ([#2981](https://github.com/cp2k/cp2k/pull/2981))
- Drop CMake option `CP2K_BUILD_DBCSR` ([#3044](https://github.com/cp2k/cp2k/pull/3044))

### Fixes

- Correct LnPP2 Basis sets. Some typos corrected Jun-Bo Lu
  ([#3100](https://github.com/cp2k/cp2k/pull/3100), [#3102](https://github.com/cp2k/cp2k/pull/3102))
- Fix Wannier localization when using LOW SPIN ROKS
  ([#3108](https://github.com/cp2k/cp2k/pull/3108))
- XAS_TDP: Fix bug leading to crashes when n_ranks>>>nex_atom
  ([#2908](https://github.com/cp2k/cp2k/pull/2908))
- EMD + Time Dependent Electric Field ([#3016](https://github.com/cp2k/cp2k/pull/3016))
- Fix parsing atom sites from CIF files ([#3092](https://github.com/cp2k/cp2k/pull/3092))
- Fix parsing of CHARMM General Force Fields ([#2956](https://github.com/cp2k/cp2k/pull/2956))

______________________________________________________________________

## 2023.2 (July 28, 2023)

- GW: Periodic open-shell and Splitting of electronic states due to spin-orbit coupling
  ([#2639](https://github.com/cp2k/cp2k/pull/2639), [#2831](https://github.com/cp2k/cp2k/pull/2831))
- GTH pseudopotential database file with spin-orbit coupling (SOC) parameters added
  ([#2848](https://github.com/cp2k/cp2k/pull/2848))
- RTP: TD Field Velocity gauge and projection TD-MOs
  ([#2623](https://github.com/cp2k/cp2k/pull/2623), [#2744](https://github.com/cp2k/cp2k/pull/2744))
- RTP: Linear density delta kick and restart ([#2543](https://github.com/cp2k/cp2k/pull/2543))
- RTP: Enabled ADMM with GAPW ([#2729](https://github.com/cp2k/cp2k/pull/2729))
- Implementation of the NVPT for APTs and AATs in velocity form
  ([#2568](https://github.com/cp2k/cp2k/pull/2568), [#2561](https://github.com/cp2k/cp2k/pull/2561))
- Intrinsic Atomic Orbitals ([#2707](https://github.com/cp2k/cp2k/pull/2707))
- Machine Learning: Add PyTorch interface, Nequip and Allegro models
  ([#2420](https://github.com/cp2k/cp2k/pull/2420), [#2528](https://github.com/cp2k/cp2k/pull/2528),
  [#2722](https://github.com/cp2k/cp2k/pull/2722))
- k-points: Implementation of the DIIS/Diag. solver
  ([#2721](https://github.com/cp2k/cp2k/pull/2721))
- TDDFPT: SOC absorption ([#2859](https://github.com/cp2k/cp2k/pull/2859))
- GAPW triplet excitation energies and forces ([#2837](https://github.com/cp2k/cp2k/pull/2837),
  [#2861](https://github.com/cp2k/cp2k/pull/2861))
- EC: Enable DC-DFT with HFX-ADMM for reference and DC calculation
  ([#2780](https://github.com/cp2k/cp2k/pull/2780))
- Add cell symmetry `HEXAGONAL_GAMMA_120` ([#2758](https://github.com/cp2k/cp2k/pull/2758))
- Grid: Rename backends, change default to `CPU` ([#2772](https://github.com/cp2k/cp2k/pull/2772),
  [#2775](https://github.com/cp2k/cp2k/pull/2775), [#2778](https://github.com/cp2k/cp2k/pull/2778))
- Grid: Enable GPU acceleration for large basis sets
  ([#2787](https://github.com/cp2k/cp2k/pull/2787), [#2793](https://github.com/cp2k/cp2k/pull/2793))
- FM: Add experimental support for [NVIDIA cuSOLVERMp](https://docs.nvidia.com/hpc-sdk/cusolvermp)
  ([#2860](https://github.com/cp2k/cp2k/pull/2860))
- Regtesting: Add `--smoketest` option ([#2501](https://github.com/cp2k/cp2k/pull/2501))
- Add support for MPI Fortran 2008 bindings ([#2486](https://github.com/cp2k/cp2k/pull/2486))
- Add support for [Apptainer/Singularity](https://apptainer.org/) containers
  ([README](https://github.com/cp2k/cp2k/blob/master/tools/apptainer/README.md))

______________________________________________________________________

## 2023.1 (January 1, 2023)

- Add gradients for SOS-MP2 and RPA incl. benchmarks
  ([#2208](https://github.com/cp2k/cp2k/issues/2208),[#2271](https://github.com/cp2k/cp2k/issues/2271),[#2473](https://github.com/cp2k/cp2k/issues/2473))
- TDDFT/Linear Response: Add GAPW/GAPW_XC and ADMM/GAPW options
  ([#2200](https://github.com/cp2k/cp2k/issues/2200))
- TDDFT: Add excited state forces as property ([#2363](https://github.com/cp2k/cp2k/issues/2363))
- RI-RPA: Allow for XC correction in ADMM RI-RPA ([#2216](https://github.com/cp2k/cp2k/issues/2216))
- RTP: Velocity gauge and magnetic delta pulse ([#2343](https://github.com/cp2k/cp2k/issues/2343))
- GW: Automatically extrapolate k-point mesh ([#2229](https://github.com/cp2k/cp2k/issues/2229))
- xTB: Add vdW options ([#2431](https://github.com/cp2k/cp2k/issues/2431))
- xTB: Fix electronic energy dependence on EPS_DEFAULT
  ([#2287](https://github.com/cp2k/cp2k/issues/2287))
- Vibrational analysis: Raman Intensities ([#2263](https://github.com/cp2k/cp2k/issues/2263))
- New pseudopotentials and basis sets ([#2472](https://github.com/cp2k/cp2k/issues/),
  [#2193](https://github.com/cp2k/cp2k/issues/2193))
- Improve NewtonX interface ([#2443](https://github.com/cp2k/cp2k/issues/2443))
- Fist: Add LAMMPS style tabulated pair potentials
  ([#2313](https://github.com/cp2k/cp2k/issues/2313))
- EC: Variational Density-Corrected DFT (DC-DFT) ([#2322](https://github.com/cp2k/cp2k/issues/2322))
- Update active space interface ([#2346](https://github.com/cp2k/cp2k/issues/2346))
- Helium: Add missing xyz output format ([#2432](https://github.com/cp2k/cp2k/issues/2432))
- SIRIUS: Add support for libvdwxc ([#2270](https://github.com/cp2k/cp2k/issues/))
- ELPA: Fix block size issue on GPU ([#2407](https://github.com/cp2k/cp2k/issues/2407))
- Drop Support for MPI 2.0 ([#2438](https://github.com/cp2k/cp2k/issues/2438))
- Add experimental CMake build system ([#2364](https://github.com/cp2k/cp2k/issues/2364))
- Fix regtests on ARM64 ([#1855](https://github.com/cp2k/cp2k/issues/1855))
- Start testing with Address Sanitizer ([#2306](https://github.com/cp2k/cp2k/issues/2306))
- Start testing on macOS Apple M1 (sponsored by [MacStadium](https://www.macstadium.com/opensource))

______________________________________________________________________

## 2022.2 (October 4, 2022)

- Minor release to fix the outdated url for Spglib in the toolchain
  ([#2262](https://github.com/cp2k/cp2k/issues/2262)).

______________________________________________________________________

## 2022.1 (July 8, 2022)

- Migrate tensor operations to new sparse matrix library DBM
  ([#1863](https://github.com/cp2k/cp2k/issues/1863))
- Add HIP support for PW ([#1864](https://github.com/cp2k/cp2k/issues/1864))
- Drop support for GCC 5 ([#1878](https://github.com/cp2k/cp2k/issues/1878))
- Add GAPW Voronoi integration ([#1919](https://github.com/cp2k/cp2k/issues/1919))
- Remove deprecated sections LIBXC and KE_LIBXC ([#1921](https://github.com/cp2k/cp2k/issues/1921))
- Add LibXC equivalents to ADMM exchange potentials
  ([#1972](https://github.com/cp2k/cp2k/issues/1972))
- Improve support for metaGGA functionals ([#1974](https://github.com/cp2k/cp2k/issues/1974))
- Use SPLA for offloading dgemm on GPUs in the mp2 module
  ([#1951](https://github.com/cp2k/cp2k/issues/1951))
- TDDFT: enable state following using transition charge finger print
  ([#1991](https://github.com/cp2k/cp2k/issues/1991))
- Add barostat for frozen atoms in absolute coordinate
  ([#2000](https://github.com/cp2k/cp2k/issues/2000))
- Fix linkage of COSMA ([#2021](https://github.com/cp2k/cp2k/issues/2021))
- Migrate to centralized `__OFFLOAD_CUDA/HIP` flags
  ([#2027](https://github.com/cp2k/cp2k/issues/2027))
- Add low-scaling SOS-Laplace MP2 forces ([#2031](https://github.com/cp2k/cp2k/issues/2031))
- Refactoring of basis set optimization code ([#2068](https://github.com/cp2k/cp2k/issues/2068))
- Add k-points for the GW self-energy ([#2073](https://github.com/cp2k/cp2k/issues/2073))
- CDFT: forces based on Hirshfeld partitioning ([#2111](https://github.com/cp2k/cp2k/issues/2111))
- RPA: Add low-scaling gradients ([#2131](https://github.com/cp2k/cp2k/issues/2131))
- MP2: Add more solvers ([#2142](https://github.com/cp2k/cp2k/issues/2142))
- GW: Add 4-center Hartree-Fock and ADMM for exchange self-energy
  ([#2145](https://github.com/cp2k/cp2k/issues/2145))
- Print vibrational modes for Newton-X ([#2146](https://github.com/cp2k/cp2k/issues/2146))
- Add partially occupied Wannier states ([#2154](https://github.com/cp2k/cp2k/issues/2154))
- Voronoi integration: Mitigated issues with symmetric structures, more diagnostic output
  ([#2171](https://github.com/cp2k/cp2k/issues/2171))
- Add GAPW_XC for TDDFPT energies ([#2178](https://github.com/cp2k/cp2k/issues/2178))

______________________________________________________________________

## 9.1 (December 31, 2021)

- Fix MacOS build ([#1316](https://github.com/cp2k/cp2k/issues/1316))
- Add NEWTONX interface ([#1794](https://github.com/cp2k/cp2k/issues/1794))
- Add Gromacs QM/MM support
  ([see also](https://manual.gromacs.org/documentation/2022-beta1/reference-manual/special/qmmm.html))
- Add experimental support for HIP and OpenCL to DBCSR
- Adopt BSD3 license for new performance critical code
  ([#1632](https://github.com/cp2k/cp2k/issues/1632))
- Add GAL21 forcefield ([#1579](https://github.com/cp2k/cp2k/issues/1579))
- Upgrade to MPI_THREAD_SERIALIZED ([#1564](https://github.com/cp2k/cp2k/issues/1564))
- Add new pseudopotentials and basis sets ([#1547](https://github.com/cp2k/cp2k/issues/1547),
  [#1551](https://github.com/cp2k/cp2k/issues/1551))
- Add analytical derivatives of the MO coefficients wrt nuclear coordinates
  ([#1706](https://github.com/cp2k/cp2k/issues/1706))
- Add forces for RI-HFX ([#1688](https://github.com/cp2k/cp2k/issues/1688))
- Add forces for TDDFT ([#1670](https://github.com/cp2k/cp2k/issues/1670),
  [#1759](https://github.com/cp2k/cp2k/issues/1759))
- Regularized RI for periodic GW ([#1776](https://github.com/cp2k/cp2k/issues/1776))
- Add beadwise constraints to PINT ([#1734](https://github.com/cp2k/cp2k/issues/1734))
- Add analytical stress tensor for NNP ([#1783](https://github.com/cp2k/cp2k/issues/1783))
- Add ghost particles and tip scan for xTB ([#1578](https://github.com/cp2k/cp2k/issues/1578))
- Add forces and stress tensor for MP2-based double-hybrids
  ([#1647](https://github.com/cp2k/cp2k/issues/1647))
- Rewrite regtesting script in Python, arguments changed slightly
  ([#1548](https://github.com/cp2k/cp2k/issues/1548))

______________________________________________________________________

## 8.2 (May 28, 2021)

- Speedup grid kernels, especially non-orthorhombic on CPU and integrate on GPU
- Upgrade to COSMA 2.5 ([#1303](https://github.com/cp2k/cp2k/issues/1303))
- Add support for ARM64
- Drop support for GCC 6 ([#1203](https://github.com/cp2k/cp2k/issues/1203))
- Upgrade to LibXC 5 and
  [harmonize its input](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL.html)
  with built in functionals
- xTB/DFTB: Add stress tensor with efield
- Motion: Add space group symmetry
- Fix multi GPU by setting the active device consistently
  ([#814](https://github.com/cp2k/cp2k/issues/814))
- Fix MOLDEN output ([#1335](https://github.com/cp2k/cp2k/issues/1335))
- XAS_TDP: Fix bug in open-shell SOC ([#1304](https://github.com/cp2k/cp2k/issues/1304))
- PINT: Fix conserved quantity in PINT-RPMD restart
  ([#1290](https://github.com/cp2k/cp2k/issues/1290))
- ELPA: As a precaution use only for large matrices by default
  ([#1444](https://github.com/cp2k/cp2k/issues/1444))
- [libvori](https://brehm-research.de/cp2k.php): Augmented .voronoi file format provides improved
  support for [TRAVIS](https://brehm-research.de/travis.php) (e.g. for spectra simulations)

______________________________________________________________________

## 8.1 (December 30, 2020)

- Fix bug affecting ADMM on GPUs ([#893](https://github.com/cp2k/cp2k/issues/893))
- Fix bug affecting Amber dihedrals ([#984](https://github.com/cp2k/cp2k/issues/984))
- Drop support for Python 2 and non-OpenMP builds
- Add interfaces to GRRM17 and SCINE codes
- Add support for Cosma (https://github.com/eth-cscs/COSMA)
- Add support for Voronoi integration of electron density (https://brehm-research.de/voronoi)
- Add support for output of electron density in compressed BQB format
  (https://brehm-research.de/bqb)
- OpenMP refactoring and speed-ups for one electron integrals
- Response code for polarizabilities: add finite difference debug, hybrid functionals, and ADMM
- Harris functional based on Kohn-Sham density
- TDDFPT code refactoring, add sTDA kernel, xTB/sTDA method
- NNP: Behler-Parrinello Neural Network Potentials
- XAS_TDP: Add OT solver and improve performance
- mGGA: Add stress tensor and fix bug ([#1116](https://github.com/cp2k/cp2k/issues/1116))
- QMMM: Add benchmarks and speedup GEEP with OpenMP
- LS: Add sign calculation based on submatrix method
- ALMO: Add trust region methods
- RI-HFX: Add resolution of identity for Hartree-Fock exchange
- RPA/GW/MP2: Several optimizations and refactoring to low-scaling implementation
- CUDA: GPU acceleration of collocate and integrate grid operations (experimental)

______________________________________________________________________

## 7.1 (December 24, 2019)

- [SIRIUS](https://github.com/electronic-structure/SIRIUS): Plane Wave module with GPU support, see
  also [this tutorial](https://www.cp2k.org/howto:running_qe_computation) for Quantum ESPRESSO
  users.
- xTB: Tight-binding module based on [](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- RPA / GW / MP2: migrated to DBCSR tensors.
- HELIUM: New canonical worm algorithm based on [](https://dx.doi.org/10.1103/PhysRevE.74.036701).
- XAS_TDP: X-ray absorption spectra simulations using linear-response TDDFT.
- NEGF: Contact-specific temperature, correct shift and scale factors.
- S-ALMO: Major refactoring, added wide variety of options.
- CDFT: Cleanup and bug fixing.
- FPGA interface for pw FFT.
- Updated libraries: DBCSR, ELPA, libint, libxc, libxsmm.
- The cp2k_shell was integrated into the main binary, simply call cp2k with `-s` or `--shell`.
- Development moved from SVN to Git

______________________________________________________________________

## 6.1 (June 11, 2018)

- Projection-operator adiabatization (POD) method
- CP2K can now do Plane Wave calculations using CPU and GPU, based on an electronic structure
  library [SIRIUS](https://github.com/electronic-structure/SIRIUS)
- Include NVIDIA P100 kernels for DBCSR
- Update toolchain
- Prevent ELPA diagonalization crashes with small matrices and/or large core counts
- Faster routines for reading and writing cube files using MPI I/O
- Docker based tests

______________________________________________________________________

## 5.1 (October 24, 2017)

- Sparse tensor framework based on DBCSR
- Flags `__ELPA2` and `__ELPA3` removed, instead use `-D__ELPA=YYYYMM` to specify library version.
- Constrained DFT, see [](#CP2K_INPUT.FORCE_EVAL.DFT.QS.CDFT) and
  [](#CP2K_INPUT.FORCE_EVAL.MIXED.MIXED_CDFT)
- Cubic scaling GW
- GW + image charge to compute electronic levels of a molecule on a metal surface
- Constraint cell optimization [](#CP2K_INPUT.MOTION.CELL_OPT.CONSTRAINT)

______________________________________________________________________

## 4.1 (October 5, 2016)

- Maximum Overlap Method (MOM)
- Modified Atomic Orbitals (MAO) Analysis
- Easier installation with an improved
  [toolchain](https://github.com/cp2k/cp2k/tree/master/tools/toolchain)
- Improved Development Tools: [prettifier](./development/code-formatting),
  [API documentation](https://apidoc.cp2k.org)
- Improved [Coding Standards](./development/coding-conventions)
- More collective variables
- [Transport with Omen:](https://www.cp2k.org/howto:cp2k_omen) improvements
- [libcp2k.h](https://github.com/cp2k/cp2k/blob/master/src/start/libcp2k.h) interface (C/C++ header)
- Remote Memory Access (RMA) for future architectures
- Various performance improvements and bug fixes
- [GTH-PBE pseudopotentials for the Lanthanide elements](https://htmlpreview.github.io/?https://github.com/cp2k/cp2k-data/blob/master/potentials/Goedecker/index.html)
- Polarized Atomic Orbitals from Machine Learning (PAO-ML)
- Cubic-scaling RPA
- Fast method for periodic ERI reducing overhead of image charge correction in QM/MM and used in
  cubic-scaling RPA
- Drop support for PLUMED 1.3 (PLUMED 2.x is now required)
- Improved linear scaling (LS) DFT MD with curvy steps
- Support for Hybrid density functionals in TDDFT

______________________________________________________________________

## 3.0 (December 22, 2015)

- Improvement of the Path integral code and use of PIGLET thermostat
- Workaround for ifort 'feature' leading to incorrectly ignored 1-4 interactions in classical MD
  simulations.
- Gradients for MP2 in the unrestricted case
- Current output for EMD
- constant E/D simulations
- RMA based DBCSR
- Improved portability (xlf90)
- G0W0 and eigenvalue self-consistent GW
- Improved testing: the make target 'test' will now regtest the code
- Basic k-point functionality for GGA DFT
- Faster TRS4 for semi-empirical runs.
- Interface to the [PEXSI library](https://pexsi.org)
- [PLUMED](https://www.plumed.org) 2.0 interface
- Interface to ELPA2015 [ELPA](https://elpa.rzg.mpg.de/)
- Filtered Basis method
- More optimized CUDA kernels
- Coupling with the quantum transport code [OMEN](https://nano-tcad.ee.ethz.ch)
- Implicit Poisson solver, with dielectric, and different boundary conditions
- Rho mixing for LS SCF enabled
- Speedup for LRIGPW method
- Support for the [ASE](https://wiki.fysik.dtu.dk/ase/) Python toolkit.
- Updated [toolchain](https://github.com/cp2k/cp2k/tree/master/tools/toolchain)
- Saveguard against known issue with CUDA cufft 7.0
- Polarized atomic orbitals
- REPEAT variant of the RESP atomic charge fitting method
- Streamlined error handling
- Support for Intel's [libxsmm](https://github.com/hfp/libxsmm/)
- Various bug fixes

______________________________________________________________________

## 2.6 (December 22, 2014)

- Allowing GPU acceleration for full matrix multiplies by using the DBCSR multiply routines
- RTP and EMD with Hartree-Fock exchange and ADMM NONE
- Improved FULL_SINGLE_INVERSE preconditioner and linear scaling PRECOND_SOLVER INVERSE_UPDATE for
  large systems
- Self-consistent continuum solvation (SCCS) model
- K-points (partial implementation for some methods)
- Improved linear scaling routines.
- Improved RPA frequency integration methods.
- QUIP Manybody potential
- Optimization of LRI basis sets
- Full Fortran 2003 compliance
- Various bug-fixes, refactoring, and speed-ups
- Collection of production grade parameters (basis-sets, pseudo-potentials, etc.)
- File discovery mechanism (controllable via compile flag `-D__DATA_DIR` or environment variable
  `$CP2K_DATA_DIR` )
- New build-system, replaced makedepf90 with novel makedep.py
- Started splitting code into sub-directories (packages) with well defined dependencies among each
  other
- Auto-tunning framework for DBCSR cuda-kernels and many readily optimized kernel-parameters
- QM/MM for DFTB
- LRIGPW
- Hirshfeld population analysis
- DM and Charge Constraint Projection based ADMM

______________________________________________________________________

## 2.5 (February 26, 2014)

- MP2 gradients and stress
- emacs plugin for input syntax highlighting (→ [CP2K Tools](https://www.cp2k.org/tools))
- vim plugin for input syntax highlighting (→ [CP2K Tools](https://www.cp2k.org/tools))
- Post-SCF linear response, including Raman
- Energy use framework for Cray
- RI-MP2 auxiliary basis optimization
- Global optimization of geometries
- SCP tight binding
- Relativistic corrections for atomic blocks
- CUDA enabled DBCSR
- Improved OMP parallelism
- Tree Monte Carlo: additional parallelism in MC
- ALMO: linear scaling for molecular systems
- Removed internal ELPA, using it as an external library instead
- Integrated molecular basis set optimization
- Langevin dynamics regions
- DCD dump option for an aligned cell (allows a reconstruction of scaled coordinates)
- and many bug fixes ...

______________________________________________________________________

## 2.4 (June 13, 2013)

- GPW-MP2 & RPA
- Adaptive QM/MM
- Non-local vDW functionals, PBEsol
- Integrated basis set optimisation
- Support for non-linear core corrected pseudos
- Additional linear scaling algorithms and properties
- Possibility of using image charges
- Periodic RESP charges
- PLUMED support
- ELPA eigensolver support
- libxc support
- Improved ifort/MKL support
- Process topology mapping for Cray Gemini

## 2.3 (Sept 03, 2012)

## 2.2 (Oct 23, 2011)

## 2.1 (Oct 6, 2010)

## 2.0 (Sep 8, 2009)
