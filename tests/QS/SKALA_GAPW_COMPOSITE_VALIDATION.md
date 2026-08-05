# SKALA GAPW composite-density validation

## Objective

Evaluate SKALA once on a GAPW density representation in which the smooth and one-center fields are
combined before SKALA constructs its nonlinear features. The target primitive fields are

```text
rho_sigma = rho_smooth_sigma
          + sum_A (rho_hard_A_sigma - rho_soft_A_sigma)

grad_rho_sigma = grad_rho_smooth_sigma
               + sum_A (grad_rho_hard_A_sigma - grad_rho_soft_A_sigma)

tau_sigma = tau_smooth_sigma
          + sum_A (tau_hard_A_sigma - tau_soft_A_sigma)
```

The model must receive these combined primitive fields before logarithms, squared-gradient
invariants, neural layers, or non-local message passing are evaluated. The existing expression

```text
E[rho_smooth] + sum_A (E[rho_hard_A] - E[rho_soft_A])
```

remains a separately testable approximation and is not the composite-density reference.

## Gradient cross terms

SKALA 1.1 constructs seven scalar point features from the two spin densities: two densities, two
squared spin-density gradients, two kinetic-energy densities, and the squared total-density
gradient. This is the expression in `src/skala/functional/model.py::_prepare_features_raw` at the
current Microsoft SKALA main commit `abc939df9a2b815e69678751a00794bdbb09309c`. For one one-center
contribution, the last feature contains the six vectors

```text
grad_rho_smooth_alpha, grad_rho_hard_alpha, -grad_rho_soft_alpha,
grad_rho_smooth_beta,  grad_rho_hard_beta,  -grad_rho_soft_beta
```

Expanding the square of their sum gives six diagonal terms and `binomial(6, 2) = 15` pairwise cross
terms. The implementation must not enumerate these products. It must sum the primitive gradient
vectors first and form the SKALA invariants afterwards, which generates all cross terms
automatically. With `s`, `h`, and `p` denoting smooth, hard, and soft and with the soft vectors
entered with a minus sign, the 15 terms are

| Pair               | Contribution                       |
| ------------------ | ---------------------------------- |
| `s_alpha, h_alpha` | `+2 grad(s_alpha) . grad(h_alpha)` |
| `s_alpha, p_alpha` | `-2 grad(s_alpha) . grad(p_alpha)` |
| `h_alpha, p_alpha` | `-2 grad(h_alpha) . grad(p_alpha)` |
| `s_beta, h_beta`   | `+2 grad(s_beta) . grad(h_beta)`   |
| `s_beta, p_beta`   | `-2 grad(s_beta) . grad(p_beta)`   |
| `h_beta, p_beta`   | `-2 grad(h_beta) . grad(p_beta)`   |
| `s_alpha, s_beta`  | `+2 grad(s_alpha) . grad(s_beta)`  |
| `s_alpha, h_beta`  | `+2 grad(s_alpha) . grad(h_beta)`  |
| `s_alpha, p_beta`  | `-2 grad(s_alpha) . grad(p_beta)`  |
| `h_alpha, s_beta`  | `+2 grad(h_alpha) . grad(s_beta)`  |
| `h_alpha, h_beta`  | `+2 grad(h_alpha) . grad(h_beta)`  |
| `h_alpha, p_beta`  | `-2 grad(h_alpha) . grad(p_beta)`  |
| `p_alpha, s_beta`  | `-2 grad(p_alpha) . grad(s_beta)`  |
| `p_alpha, h_beta`  | `-2 grad(p_alpha) . grad(h_beta)`  |
| `p_alpha, p_beta`  | `+2 grad(p_alpha) . grad(p_beta)`  |

The first six already occur in the two spin-resolved squared-gradient features. The remaining nine
couple alpha and beta and additionally occur in the squared total-density gradient. Thus `15` counts
all distinct component pairs in the complete composite field, whereas `9` counts only the cross-spin
subset. Counting only the pairs newly introduced relative to the smooth baseline gives `14`:
`s_alpha, s_beta` is already present in `|grad(s_alpha + s_beta)|^2`. Thus `15` denotes all distinct
component pairs in the complete field, whereas `14` denotes the additional pairs introduced by the
hard-minus-soft reconstruction. The implementation is independent of this counting convention
because it combines the primitive vectors before constructing the invariant.

## Backend scope

| Backend and representation                                     | Composite correction required                        |
| -------------------------------------------------------------- | ---------------------------------------------------- |
| Molecular GauXC, all-electron `METHOD GAPW`                    | No; GauXC evaluates the original AO density directly |
| Molecular GauXC, GTH/ECP direct valence density                | No                                                   |
| Molecular GauXC followed by a GAPW one-center SKALA correction | Yes                                                  |
| Native-grid GPW or GAPW kinds treated as GPW                   | No                                                   |
| Native-grid all-electron GAPW                                  | Yes                                                  |
| Native-grid pseudopotential GAPW one-center representation     | Yes                                                  |
| Native-grid `METHOD GAPW_XC` with one-center terms             | Yes                                                  |

For molecular GauXC, CP2K currently creates the GauXC basis from the original `ORB` basis. Any
additional hard-minus-soft SKALA correction must therefore be justified against that already-direct
AO density evaluation rather than assumed to be a smooth-density contribution.

For pseudopotential GAPW, the public representation selector has four deliberately distinct values:

| Value                  | Meaning                                                                      |
| ---------------------- | ---------------------------------------------------------------------------- |
| `DIRECT_VALENCE`       | Evaluate SKALA on the direct valence density; this is the default.           |
| `PAW_ONE_CENTER`       | Combine smooth, hard, and soft primitive fields before SKALA evaluation.     |
| `PAW_ONE_CENTER_SPLIT` | Retain the separately evaluated legacy expression as an explicit diagnostic. |
| `CP2K_DEFAULT`         | Recover the representation selected by the pre-existing CP2K kind settings.  |

`PAW_ONE_CENTER` includes the complete nonlinear coupling by construction. It is not an alias for
the split expression, and it does not enumerate or invoke the model separately for the 15 gradient
component-pair terms.

## Prototype stages

1. Build a common-grid reference that combines smooth, hard, and soft primitive fields before one
   SKALA forward call. Start with an energy-only CPU path and retain the existing split expression
   for side-by-side diagnostics.
1. Compare the combined native result with molecular GauXC all-electron references for H2, H2O, and
   Ar. Use identical geometries, orbital bases, SCF thresholds, and converged quadratures.
1. Validate the returned density, gradient, and kinetic-energy-density derivatives against finite
   differences of the combined energy.
1. Validate nuclear forces and the finite-system virial by finite differences before enabling an
   analytical periodic stress path.
1. Batch independent one-center work, reuse tensor buffers, and distribute the final data flow over
   MPI ranks and GPUs only after the combined expression is validated.

## Common FFT-grid prototype evidence

All values below use the same converged molecular GauXC density matrix in each paired comparison.
Atoms are centered inside a non-periodic cell; earlier geometries with atoms on or outside the cell
boundary are not valid reference data.

For H2 at 400 Ry, increasing the cubic cell reduces the native-minus-GauXC XC error from
`0.672 kJ mol^-1` at 6 bohr through `0.485 kJ mol^-1` at 8 bohr to `0.203 kJ mol^-1` at 10 bohr. The
GauXC XC energy is invariant at fixed density matrix. The native-grid quadrature depends directly on
the finite cell and grid, whereas the CP2K analytic Poisson solver affects the Hartree and total
energies but not the GauXC XC quadrature. At a fixed 6-bohr cell, raising the cutoff from 400
through 800 to 1200 Ry reduces the native XC error from `0.672` through `0.626` to
`0.558 kJ mol^-1`.

The same global-grid reconstruction is not viable for heavier all-electron atoms at practical
cutoffs. The fixed-density XC energies are:

| System | Molecular GauXC (Ha) | Composite 200 Ry (Ha) |   400 Ry (Ha) |    800 Ry (Ha) |  1200 Ry (Ha) |
| ------ | -------------------: | --------------------: | ------------: | -------------: | ------------: |
| H2O    |        -9.2654879381 |         -5.0158042290 | -8.5573136070 |  -9.3209592495 | -9.2801852100 |
| Ar     |       -30.8749302599 |        -28.6311544280 |       VXC NaN | -30.7235135017 |       VXC NaN |

An independent current-stack repeat from a newly converged molecular GauXC H2O wavefunction gives
`-9.2657876989 Ha` for molecular GauXC and `-9.2065996351`, `-9.3211507221`, and `-9.2804356591 Ha`
for the common-grid reconstruction at 400, 800, and 1200 Ry. The corresponding errors are `+155.40`,
`-145.36`, and `-38.46 kJ mol^-1`; the sequence is strongly oscillatory rather than monotonically
convergent. The reconstructed electron and kinetic-energy-density integrals stay fixed at
`10.000000006` and `75.614911724 Ha`, but the 1200-Ry fields still reach negative minima of
`-8.14e-3` for rho and `-3.55 Ha` for tau. Direct AO collocation at 800 Ry is likewise unresolved:
it integrates to `10.130304` electrons and `69.028741 Ha` of kinetic-energy density and gives
`-10.2540475089 Ha`. The common-grid 1200-Ry model batch also exceeds a 46-GB A40 after using about
`42 GB` and requesting another `7.2 GB`. These current-stack controls reinforce that a global FFT
grid is a diagnostic limit, not the production all-electron representation.

For Ar at 400 Ry the reconstructed fields have `rho_min = -7.344`, `rho_max = 143.524`,
`tau_min = -3614.694`, and `tau_max = 21331.536`. The exact electron and kinetic-energy integrals
therefore do not imply a pointwise physical representation: Fourier truncation produces severe
ringing. Direct full-AO collocation on the same 200-Ry grid integrates to 40.529 electrons for Ar
instead of 18, which independently exposes the unresolved core length scale. Clipping these fields
would change the functional and is not an acceptable remedy.

The existing split expression retains atom-centered core resolution but remains non-identical for
SKALA. Adding the regular-grid and printed one-center XC terms gives errors relative to molecular
GauXC of approximately `1.00 kJ mol^-1` for H2, `8.74 kJ mol^-1` for H2O, and `-2.72 kJ mol^-1` for
Ar at 200 Ry. Its occasionally smaller error than the common FFT-grid result does not make the split
expression exact.

For the one-atom Ar limit, the hard one-center term can be compared directly because no interatomic
or overlapping-augmentation contribution exists. A `SUPERFINE/UNPRUNED` GauXC grid gives
`-30.8745824855 Ha`. The CP2K hard one-center grid gives `-30.8745832715 Ha` with 120 radial and at
least 120 requested Lebedev points, a difference of only `-0.0021 kJ mol^-1`. At 150/150 the
difference is `0.0091 kJ mol^-1`. The earlier 50/50 result differs by several kJ/mol; its angular
and radial convergence tests show that the radial grid is the limiting component. Thus the
atom-centered primitive representation can reproduce molecular GauXC AE accuracy without an extreme
global plane-wave cutoff.

The production path must therefore preserve atom-centered core resolution while presenting a single
combined primitive-field data set to SKALA. A suitable hybrid quadrature has regular-grid rows in
the smooth/interstitial region and atom-centered rows in core regions. On every retained row it must
form `smooth + hard - soft` for density, gradient, and kinetic-energy density before feature
construction. Every atom block must contain all of its rows and complete primitive fields so that
local gradient cross terms and the atom-block non-local layers are retained. Complete atom blocks
may be evaluated independently because Skala 1.1 does not reduce across its atom dimension. The
quadrature partition must be differentiable with respect to atom positions and strain before
analytical forces or virials can be enabled.

## Atom-centered composite prototype evidence

The second prototype interpolates the smooth native-grid fields to the existing GAPW radial/Lebedev
grids, adds the one-center hard-minus-soft fields there, and applies a smooth atom partition. The
initial reference sent every atom's rows through one SKALA call; the distributed implementation
evaluates complete atom blocks independently with equivalent model semantics. Unlike the common
FFT-grid prototype, this preserves core resolution. All comparisons use the same restart density
matrix and the same geometry on both sides.

| System | Molecular GauXC (Ha) |    200 Ry (Ha) |    400 Ry (Ha) |    800 Ry (Ha) |   1200 Ry (Ha) |
| ------ | -------------------: | -------------: | -------------: | -------------: | -------------: |
| H2     |        -0.6846515268 |  -0.6847252974 |  -0.6845939964 |  -0.6846837102 |  -0.6846569988 |
| H2O    |        -9.2654879804 |  -9.2645844756 |  -9.2641194480 |  -9.2636593880 |  -9.2636472827 |
| Ar     |       -30.8745824855 | -30.8765301041 | -30.8749614247 | -30.8747078359 | -30.8746458060 |

At 1200 Ry, H2 and Ar differ from their molecular references by `0.014` and `0.166 kJ mol^-1`,
respectively. The Ar electron integral is `18.0000054`. Linear native-grid interpolation converged
far more slowly and still missed the Ar reference by about `8 kJ mol^-1` at 1200 Ry. These energy
tests established cubic interpolation as the minimum useful order; the active force-consistent path
now uses the higher-order operator described below.

H2O exposes a separate reconstruction limit. At 800 Ry the 50/50 and 120/120 radial/Lebedev grids
give `-9.2636593880` and `-9.2636036010 Ha`, and increasing `GAPW_1C_BASIS` from `EXT_SMALL` to
`EXT_VERY_LARGE` changes the result only to `-9.2638129396 Ha`. Enlarging and recentering the cell
from 5 to 8 Angstrom gives `-9.2638384578 Ha`; vacuum truncation is therefore a secondary effect.
The remaining approximately `1.65 mHa` (`4.33 kJ mol^-1`) difference is a finite GAPW reconstruction
error rather than FFT core ringing.

On the current stack, a centered 400-Ry H2O repeat gives an atom-composite XC energy of
`-9.26397492987918 Ha`, compared with `-9.26578769888037 Ha` from molecular GauXC at the same
density, a difference of `4.76 kJ mol^-1`. The reconstructed fields integrate to `9.9998139126`
electrons and `75.5903088134 Ha` of kinetic-energy density. Their extrema are `rho_min = -9.33e-6`,
`rho_max = 146.34`, `tau_min = -8.23e-6`, and `tau_max = 4683.44`. The native XC construction takes
`2.10 s` and reaches approximately `1.47 GiB` resident memory in this test.

A fixed-density molecular GauXC control cleanly separates the XC integration from CP2K's isolated
electrostatics. Reusing the same H2O wavefunction in 5, 8, and 10 Angstrom non-periodic cells with
`POISSON_SOLVER ANALYTIC` gives the identical GauXC-SKALA XC energy `-9.26548788983717 Ha` in all
three cells. The corresponding total CP2K energies are `-76.4405994166`, `-76.4091589095`, and
`-76.4090575826 Ha`, respectively. GauXC therefore sees neither the cell nor the Poisson solver at
fixed AO density, while the Hartree/electrostatic part of the CP2K calculation remains substantially
unconverged in the 5 Angstrom cell. A native-grid XC comparison may additionally depend directly on
the cell through its finite grid and interpolation.

The component diagnostic confirms that the correction is physically active. For H2O at 800 Ry, the
smooth field alone integrates to `7.9044306` electrons and gives `-4.1371417070 Ha`. Adding only the
hard-minus-soft density and gradient restores `9.9996853` electrons and gives `-9.6318609786 Ha`;
adding the tau correction yields the full `-9.2636593880 Ha`. Explicitly adding the correction of
atom B on atom A's grid does not change these H2O values at the shown precision, because the
augmentation support is smaller than the interatomic grid overlap in this case.

The Torch backward derivative of the composite energy was checked by scaling density, gradient, and
tau together and comparing its contraction with a central finite difference at step `3e-3`:

| System | Analytic contraction (Ha) | Finite difference (Ha) | Difference (Ha) |
| ------ | ------------------------: | ---------------------: | --------------: |
| H2     |             -0.9166066436 |          -0.9166063366 |        -3.07e-7 |
| H2O    |            -12.5420315586 |         -12.5420305881 |        -9.70e-7 |
| Ar     |            -41.4126471212 |         -41.4126386975 |        -8.42e-6 |

The three primitive channels were also scaled separately. At the same finite-difference step, the
largest component error is about `1.4e-6 Ha` for H2, `2.0e-5 Ha` for H2O, and `1.7e-4 Ha` for Ar.
These larger heavy-atom values are the expected central-difference truncation error of the
individual nonlinear directions; the combined density-matrix scaling direction above remains the
physical VXC acceptance check.

The complete VXC adjoint was then projected through both CP2K representations. For the smooth field,
the transpose of the same native-grid interpolation was followed by CP2K's discrete gradient
divergence and tau potential construction. For the one-center field, the same model derivative was
applied to hard and soft channels, integrated with `gaVxcgb_GC` and `dgaVtaudgb`, and combined as
hard minus soft. The padded old one-center integral basis was mapped to the compact density-matrix
basis with the CP2K `n2oindex` transformation before contraction. The resulting adjoint errors are:

| System | Smooth PW adjoint error (Ha) | One-center matrix adjoint error (Ha) |
| ------ | ---------------------------: | -----------------------------------: |
| H2     |                    -2.29e-14 |                             1.75e-14 |
| H2O    |                    -7.73e-14 |                             1.03e-13 |
| Ar     |                     6.75e-14 |                            -6.11e-13 |

For H2O, separate density, gradient, and tau matrix contractions also agree with their primitive
field contractions to approximately `3e-14 Ha` or better. This validates the model VXC, the cubic PW
interpolation adjoint, the gradient-divergence pair, the hard-minus-soft sign, spin factors,
one-center basis indexing, and the tau projection. The model and matrix VXC derivatives are
therefore validated independently of nuclear-coordinate and strain derivatives. Those latter
derivatives are assessed below.

The atom-centered route was subsequently connected as the active XC energy and potential rather than
an energy-only diagnostic. The smooth native-grid XC energy is then set to zero, the combined energy
is stored as the GAPW one-center contribution, and the same Torch backward pass supplies both the
PW-grid and one-center VXC adjoints. With a 400-Ry grid, H2 converges at `EPS_SCF 1e-8` in one
restart step to `-1.167484249050976 Ha`; the combined-model energy `-0.6968038039306965 Ha` and
CP2K's stored GAPW XC contribution `-0.69680380393070 Ha` agree to the printed precision. H2O
converges with OT (`STEPSIZE 0.1`, `EPS_SCF 1e-7`) in nine steps to `-76.400754306475847 Ha`, and Ar
with a 120/120 radial/Lebedev grid converges in five steps to `-527.550455333526088 Ha`.

For each converged density, a one-step molecular GauXC calculation with the same wavefunction,
geometry, cell, and CP2K electrostatics isolates the residual representation error:

| System | Atom-composite XC (Ha) | Molecular GauXC XC (Ha) | Difference (kJ/mol) |
| ------ | ---------------------: | ----------------------: | ------------------: |
| H2     |    -0.6968038039306965 |       -0.69675310757959 |              -0.133 |
| H2O    |    -9.2983821178016015 |       -9.29950902690625 |               2.959 |
| Ar     |    -30.874957688168386 |      -30.87457906362269 |              -0.994 |

These differences are present in the XC contribution itself at fixed density. They are therefore
finite GAPW reconstruction/quadrature errors rather than effects of the Poisson solver, vacuum, or
SCF convergence. The active VXC reaches a stationary SCF solution for all three systems. Analytical
coordinate derivatives of the hybrid atom-centered quadrature and strain derivatives of the periodic
regular-grid representation are assessed below.

The hard-minus-soft reconstruction conserves particle number at the density-representation level;
the finite radial/Lebedev quadrature does not integrate every represented density exactly. The
coarse H2/GTH regression grid gives `1.99963990` electrons, while 480/60 Ry with 120 radial and 120
requested Lebedev points gives `1.99969429` electrons. The corresponding all-electron values include
`9.9998139` electrons for H2O and `18.0000054` for Ar. These residuals are numerical quadrature and
interpolation errors, not a change in the CP2K density-matrix trace. No density-dependent
renormalization is applied because it would define a modified functional and introduce additional
VXC, force, and virial derivatives. The H2/GTH regression checks the composite integral directly
against two electrons with an explicit `5e-4` quadrature tolerance.

## Atom-composite force evidence

The atom-composite force combines four explicit coordinate derivatives from the same Torch backward
pass: the model atom coordinates, the model grid coordinates, interpolation of the smooth PW fields
onto a moving atom grid, and the smooth atom-partition weights. The interpolation derivative is
evaluated analytically. A direct central finite difference of the partition weight agrees with its
analytic derivative to approximately `1e-11 Ha/bohr` for centered H2.

An initial H2 force error of approximately `1.36e-2 Ha/bohr` was traced to CP2K's
`GAPW_ACCURATE_XCINT` weight force. That term differentiates the Gaussian hard/soft weights of the
split GAPW XC quadrature. The atom-composite reference replaces that quadrature by partitioned
radial/Lebedev weights, so adding the split-quadrature derivative differentiates an energy that is
no longer present. `GAPW_ACCURATE_XCINT T` remains explicit in the validation inputs, but the
combined `PAW_ONE_CENTER` SKALA term is independent of its legacy hard/soft integration weights and
their force and virial derivatives. The keyword retains its normal meaning for classical GAPW XC,
`CP2K_DEFAULT`, and the explicitly split diagnostic.

With `GAPW_ACCURATE_XCINT T`, a 400-Ry grid, a `1e-3 bohr` central difference, and the molecule
centered in a non-periodic 6-Angstrom cell, the active H2 results are:

| Atom-composite fields | Numerical force (Ha/bohr) | Analytical force (Ha/bohr) | Difference (Ha/bohr) |
| --------------------- | ------------------------: | -------------------------: | -------------------: |
| Smooth only           |               -0.03706205 |                -0.03706377 |             -1.72e-6 |
| Smooth + hard - soft  |               -0.03831412 |                -0.03831689 |             -2.77e-6 |

The smooth-only result validates the model and PW interpolation force chain. The full result adds
the one-center density, gradient, and kinetic-energy-density matrix adjoints and demonstrates that
they preserve the force agreement.

The original four-point cubic interpolation nevertheless left a translational eggbox force of
approximately `7.1e-4 Ha/bohr` per H atom at 200 Ry and `4.3e-4 Ha/bohr` at 480 Ry. The active
operator therefore uses eight Lagrange points in every grid direction. Its value, Cartesian
derivative, and transpose are generated from the same weights, preserving the discrete energy/VXC
adjoint. For a centered H2/GTH `PAW_ONE_CENTER` calculation the transverse force decreases to
`7.18e-5 Ha/bohr` per atom at 200 Ry and `7.30e-6 Ha/bohr` at 480 Ry. At 480 Ry the total energy is
`-1.167840129312483 Ha`; the analytical force on atom 2 along the bond is `7.13740923e-3 Ha/bohr`,
while a `1e-4 Angstrom` central energy difference gives `7.13998555e-3 Ha/bohr`, an absolute
difference of `2.58e-6 Ha/bohr`.

At the same 480-Ry density, the complete feature-space VXC contraction and finite difference differ
by `4.13e-7 Ha`. The smooth PW interpolation adjoint agrees to `4.7e-15 Ha`, and the one-center
field/matrix contractions agree to about `4e-15 Ha`. This independently verifies that increasing the
interpolation order did not alter the energy/VXC duality.

Independent full-composite checks extend the force validation beyond H2:

| System and cutoff     | Numerical force (Ha/bohr) | Analytical force (Ha/bohr) | Absolute difference (Ha/bohr) |
| --------------------- | ------------------------: | -------------------------: | ----------------------------: |
| H2, 400 Ry            |               -0.03831366 |                -0.03831689 |                       3.22e-6 |
| H2O, 200 Ry           |               -0.00308539 |                -0.00313175 |                       4.64e-5 |
| Ar, 200 Ry            |               -0.00568936 |                -0.00551046 |                       1.79e-4 |
| Ar, 400 Ry            |                0.00059722 |                 0.00057820 |                       1.90e-5 |
| compressed H2, 200 Ry |                2.41426855 |                 2.41461384 |                       3.45e-4 |

The Ar absolute error decreases by almost one order of magnitude from 200 to 400 Ry. Its residual
translation force changes sign as the grid is refined, which identifies finite-grid translational
invariance rather than a missing model derivative as the limiting error. The compressed H2 check
forces the one-center augmentation regions to overlap and still remains within the diagnostic
`5e-4 Ha/bohr` tolerance.

## Virial and boundary-condition evidence

A cell finite difference is not a valid molecular-virial reference for a non-periodic CP2K
calculation with `POISSON_SOLVER ANALYTIC`. The same debug calculation already fails for ordinary
PBE: at 200 Ry its numerical diagonal virial is `(-0.00411853, -0.00411850, -0.04344512) Ha`, while
the analytical total virial is `(-0.29239356, -0.29239356, -0.53228850) Ha`. This discrepancy is
independent of SKALA and must not be used to judge the atom-composite molecular virial.

The appropriate isolated-system test scales the nuclear displacements about the molecular center
while keeping the CP2K cell fixed. For H2 at 400 Ry, the total-energy finite difference gives
`dE/ds = 0.0579297422 Ha`; contraction of the analytical forces gives `0.0579308974 Ha`, a
difference of `1.16e-6 Ha`. Thus the molecular atom-composite force virial is validated without a
cell derivative of the analytic 0D Poisson solver.

For a fully periodic 6-Angstrom H2 cell, the corresponding PBE control agrees: the numerical and
analytical diagonal virials differ in their summed components by `3.52e-6 Ha` at 200 Ry. The
periodic regular-grid composite SKALA path gives an analytical diagonal virial of approximately
`(-0.00055585, -0.00055585, -0.05822297) Ha`. With a `1e-3 bohr` absolute cell displacement, the y
and z components agree with finite differences to `4.0e-6` and `5.5e-6 Ha`; the x finite difference
is dominated by an approximately `9e-8 Ha` energy noise because the expected pair-energy change is
only about `1e-7 Ha`. Increasing the displacement to `5e-3 bohr` makes the x component agree to
`2.8e-9 Ha`, while y and z then incur finite-step/grid errors of approximately `1.2e-4 Ha`.

A full nine-component check in a skewed periodic cell initially appeared to expose a native-SKALA
strain error. At 200 Ry and a `0.003 bohr` cell-matrix displacement, pure GPW/SKALA gave a maximum
`PV` difference of `1.276e-3 Ha` and a summed absolute difference of `1.038e-2 Ha`. Increasing the
step to `0.03 bohr` reduced these values almost inversely with the step, to `1.296e-4` and
`1.049e-3 Ha`. Converting every stored model parameter and buffer to double precision and replacing
the traced float casts did not change the direct `h11`, `h12`, or `h22` derivative errors at the
small step. The effect was therefore neither a missing analytical term nor float32 model arithmetic.

Layout diagnostics identified the discontinuity. The smooth native-grid partition omits atom rows
whose partition weight is at most `1e-12`, but each retained row previously carried the full native
grid weight in SKALA's internal non-local atom integration. Under the `h11` displacement, the row
count changed from `313530` for the plus deformation to `313540` for the minus deformation. The
corresponding `h12` counts were `313538` and `313534`. Rows crossing a nominally negligible outer
partition threshold therefore entered or left the non-local moments with finite weight.

The corrected layout retains the sparse cutoff but multiplies the internal atom-grid weight by a
quintic taper. The taper is zero at the `1e-12` inclusion boundary, reaches one at `1e-11`, and has
zero first derivative at both endpoints. Its derivative with respect to the smooth partition weight
is included in both atom forces and strain derivatives through the Torch derivative with respect to
`atomic_grid_weights`. Rows may still enter or leave the sparse layout, but their energy and first
derivative vanish continuously at that boundary. The unstrained H2 energy changes by only
`5.0e-6 Ha` (`0.013 kJ mol^-1`).

With the taper, all nine direct cell-matrix derivatives at 200 Ry and `0.003 bohr` agree within
`5.80e-6 Ha/bohr`. The maximum and summed absolute `PV` differences are `6.42e-5` and `1.64e-4 Ha`,
respectively, compared with `1.50e-4` and `3.29e-4 Ha` for an ordinary TPSS control at the same
cutoff and step. Raising the TPSS control to 400 Ry reduces its maximum difference to `7.76e-6 Ha`,
confirming the expected finite-grid convergence of CP2K's common density-gradient and
kinetic-energy-density virial machinery. Changing `CELL_REF` together with the deformed cell alters
the tested SKALA energies by less than `6e-10 Ha`; it is not the source of the derivative residual.
The direct `h11` error decreases from `5.80e-6` through `2.82e-6` to `1.27e-6 Ha/bohr` when the cell
step is increased from `0.003` through `0.006` to `0.012 bohr`. No step-independent missing virial
term is visible.

The same full nine-component test was repeated for the all-electron GAPW common-grid composite. At
200 Ry and a `0.003 bohr` cell-matrix displacement, every direct derivative agrees within
`6.24e-6 Ha/bohr`; the maximum and summed absolute `PV` differences are `6.56e-5` and `2.66e-4 Ha`,
respectively. The unstrained total energy is `-1.147863233501682 Ha`. This verifies that the
smooth-partition and taper strain chain applies unchanged when the native smooth fields are combined
with the all-electron hard-minus-soft reconstruction before SKALA feature formation.

The same taper derivative enters nuclear forces. In a separate periodic H2 coordinate-scaling test
at 200 Ry, force contraction gives `dE/ds = 0.1006325638043 Ha`, while a `1e-3` central energy
difference gives `0.1006324635400 Ha`. Their difference is `1.00e-7 Ha`. This independently
validates the derivative of the tapered internal atom-grid weight with respect to nuclear
coordinates.

The public `PAW_ONE_CENTER` selector was then checked for periodic pseudopotential GAPW cases with
`GAPW_ACCURATE_XCINT T`. For H2/GTH at 200 Ry, a Cartesian displacement of atom 2 gives
`F_z = -0.0398091749 Ha/bohr` from a central energy difference and `-0.0398061127 Ha/bohr`
analytically, an absolute difference of `3.06e-6 Ha/bohr`. For HCl/ccECP, the force projected onto
the skew cell's third lattice vector has the following grid convergence:

| Cutoff | Numerical force (Ha/bohr) | Analytical force (Ha/bohr) | Absolute difference (Ha/bohr) |
| ------ | ------------------------: | -------------------------: | ----------------------------: |
| 240 Ry |                0.02763438 |                 0.02764291 |                       8.53e-6 |
| 480 Ry |                0.02763932 |                 0.02763995 |                       6.28e-7 |

The ECP error decreases by more than a factor of 13 when the cutoff is doubled. The remaining 240-Ry
difference is therefore finite-grid error rather than a missing one-center force term. The
corresponding skew-cell stress checks compare direct cell-matrix derivatives. At 200 Ry the H2/GTH
`h11`, `h12`, and `h22` errors are `2.54e-6`, `2.14e-7`, and `5.23e-7 Ha/bohr`, respectively. At 240
Ry the HCl/ccECP `h11` error is `5.92e-8 Ha/bohr`.

## NLCC evidence

NLCC enters the combined representation before Skala feature construction. CP2K augments the
primitive density and density gradient, while the kinetic-energy density remains valence-only. In
the molecular `PAW_ONE_CENTER` path, GTH or SGP core Gaussians and their Cartesian first and second
derivatives are evaluated directly on the partitioned atom-centered quadrature. The same Skala
feature adjoints then provide both the moving-grid and pseudopotential-center force terms. This
avoids interpolating the sharply localized core density from a finite PW grid and does not apply a
density rescaling.

For molecular HF with GTH-NLCC-PBE potentials, `GAPW_ACCURATE_XCINT T`, a 400-Ry grid, and a
`1e-3 bohr` displacement of the F atom along the bond, the numerical and analytical forces are
`-0.02781736` and `-0.02779731 Ha/bohr`, respectively. Their absolute difference is
`2.01e-5 Ha/bohr`. The residual total translation force decreases from `3.57e-3 Ha/bohr` at 400 Ry
to `1.88e-3 Ha/bohr` at 800 Ry. A no-NLCC control shows the same transverse finite-grid behavior, so
it is not evidence of a missing NLCC derivative.

For periodic native-grid HF/GTH-NLCC at 400 Ry, the numerical diagonal virial components are
`(0.04050825, 0.02148169, 0.04799977) Ha`, compared with analytical values of
`(0.04054365, 0.02157227, 0.04785084) Ha`. The component errors are
`(-3.54e-5, -9.06e-5, 1.49e-4) Ha`, or at most `0.42%`. A separate coordinate finite difference
gives an absolute force difference of `3.65e-4 Ha/bohr` at the same cutoff. These residuals are
within the finite-grid diagnostic tolerance and converge from the substantially coarser 200-Ry case.

The atom-centered molecular composite cannot be applied unchanged to a periodic cell: its
radial/Lebedev grids extend over periodic images, and the present nearest-image smooth partition
would count the periodic density repeatedly. A periodic atom-centered implementation requires an
explicit lattice-image partition. Until that exists, periodic composite evaluation uses the common
regular grid, while the atom-centered route remains molecular.

The boundary-condition control is separate from both derivative tests. At fixed AO density matrix,
GauXC receives no CP2K cell or Poisson-solver information and its XC energy is exactly invariant to
vacuum. Cell-size dependence of a self-consistent GauXC calculation can arise only indirectly via
CP2K's Hartree/electrostatic potential and the resulting density matrix. Native-grid XC additionally
depends directly on the finite grid and interpolation.

## MPI reconstruction and derivative reductions

The common-grid reconstruction now receives the density structure selected by the caller explicitly.
This distinction is essential for `METHOD GAPW_XC`: CP2K builds its smooth density and
kinetic-energy density in `rho_xc`, whereas the ordinary `rho` structure does not contain a valid
tau field in that branch. Reading `qs_env%rho` implicitly produced a rank-dependent soft tau and
placed the GTH test energy near `-0.4188 Ha`. Passing `rho_xc` through the reconstruction restores
the expected direct-valence region, `-0.9656916972 Ha` in the two-rank smoke test. The analogous
HCl/ccECP value is `-15.3605678595 Ha`.

The discrete one-center adjoint also keeps translation forces and strain derivatives rank-local.
Quickstep globally reduces `rho_elec` and every virial component after assembling all force terms;
reducing these contributions inside the reconstruction as well counted them once per MPI rank. Only
the one-center potential matrices are reduced immediately, because the subsequent atom-distributed
`update_ks_atom` path requires their global values. For a fully converged, deliberately coarse 80-Ry
GAPW_XC/GTH calculation, the one- and two-rank energies differ by `5.9e-9 Ha`; the bond- force
difference is `1.75e-5 Ha/bohr`, and the largest diagonal stress difference is about `9 bar`. Before
removal of the duplicate reduction, the corresponding force and stress differences were
approximately `7.4e-3 Ha/bohr` and `1.8e4 bar`.

The molecular atom-grid composite is also collective. Target atoms and their cross-atom corrections
are assembled only on their owning ranks. The derived radial hard/soft fields are packed and
replicated once per kind because a local target grid can overlap the augmentation field of an atom
owned elsewhere. Each rank then evaluates its complete local atom blocks and the communicator sums
the rank-local energies. Density, gradient, tau, coordinate, and weight adjoints remain on the
owning rank for local one-center matrix projection and force assembly. The smooth PW interpolation
adjoint is the remaining full-grid object and is summed before its transpose projection.

This atom decomposition is exact for the official Skala 1.1 architecture. In
`src/skala/functional/model.py`, every tensor before the final energy sum retains an independent
atom dimension. The so-called non-local layers contract fine-grid rows only into the corresponding
atom's coarse feature; they do not reduce across atoms. Microsoft's production PySCF integration
therefore also evaluates complete atom chunks independently in `src/skala/pyscf/numint.py` and sums
their energies and VXC. CP2K follows the same decomposition. It is crucial that each local atom
block contains the full smooth plus hard-minus-soft primitive density, including overlapping
augmentation tails, before this model evaluation.

A direct three-atom TorchScript test independently verifies this model property for both official
Skala 1.1 checkpoints. Splitting the complete input into a one-atom and a two-atom block leaves the
summed energy unchanged to the last printed bit. After backward propagation, the largest absolute
difference over the adjoints of density, density gradient, kinetic-energy density, fine-grid
coordinates and weights, atomic-grid weights, and coarse atom coordinates is `5.55e-17`. The test
loads the CUDA checkpoint on the CPU only after applying the same device-constant remapping as the
CP2K Torch adapter; it therefore validates the serialized model graph rather than CUDA arithmetic.

With the root-evaluated reference implementation, the compact H2 `ENERGY_FORCE` regression gives
exactly `-0.967804782209128 Ha` with one, two, and four MPI ranks. After switching to rank-local
model calls, a macOS diagnostic with the same official checkpoint gives `-0.967804785614801 Ha` with
one rank and `-0.967804783717418 Ha` with two and four ranks. The `1.90e-9 Ha` difference is the
floating-point summation order of two independent atom energies. The largest Cartesian force change
is below `6e-10 Ha/bohr`; two and four ranks are identical, including ranks with empty atom blocks.

The final Linux run on Terok confirms the same behavior with two OpenMP threads per rank. The one-,
two-, and four-rank energies are `-0.967804782209128`, `-0.967804782473838`, and
`-0.967804782473838 Ha`, respectively. Their range is `2.65e-10 Ha`, and the largest printed force
component variation is approximately `2e-9 Ha/bohr`. The corresponding wall times are `14.83`,
`14.25`, and `15.59 s`; this two-atom case is dominated by model and process startup and is not a
scaling benchmark.

An independent AArch64 CPU build on Spark, linked against CPU-only DBCSR, gives `-0.967804784260037`
and `-0.967804784262238 Ha` with one and two MPI ranks. The energy range is `2.20e-12 Ha`, and the
largest printed force-component variation is approximately `1.1e-9 Ha/bohr`. Relative to the Terok
result, the energy and force differences remain below `2.1e-9 Ha` and `7e-9 Ha/bohr`, respectively.
This check uses one OpenMP thread per MPI rank: the available AArch64 installation combines separate
LibTorch and host BLAS/OpenMP runtimes and is not stable with nested two-thread execution. That
installation constraint is independent of the rank-local composite-density result.

An H2O/GTH diagnostic exercises unequal atom blocks and two kinds. Its one-, two-, three-, and
four-rank energies are `-16.485390778331052`, `-16.485390777122827`, `-16.485390778588240`, and
`-16.485390778588240 Ha`, respectively. The complete decomposition range is `1.47e-9 Ha`; the
largest printed Cartesian-force variation is about `6e-9 Ha/bohr`. Three and four ranks are
identical because the fourth rank owns no atom block.

The periodic HCl/ccECP stress smoke test covers unequal one-center grids and the ECP adjoint path.
Its one- and two-rank energies are `-13.373000499411184` and `-13.373000499411177 Ha`; the largest
printed Cartesian-force component changes by about `1.8e-8 Ha/bohr`. The analytical stress traces
differ by approximately `5 bar` against an absolute magnitude of `4.57e6 bar` in this deliberately
coarse diagnostic.

The final validation is rebased on CP2K master commit `50ddb19844cafd1de3e2dc9d76467b0f900461b2` and
uses GauXC commit `ea0f1fc7c02497aa950376eda61686ba30256999`. This master revision routes native
SKALA atom chunks by default for energy and VXC while retaining the full differentiable graph for
force and stress evaluations. After removing a redundant regular-grid SKALA evaluation from the
atom-composite route, the final production-driver results on Terok are:

| Directory                 | Layout                         | Matcher result | Driver wall time | Slow tests |
| ------------------------- | ------------------------------ | -------------: | ---------------: | ---------: |
| `regtest-gauxc-gapw-gth`  | 2 MPI x 2 OpenMP, CPU          |          13/13 |          51.01 s |          0 |
| `regtest-gauxc-gapw-ecp`  | 2 MPI x 2 OpenMP, CPU          |            8/8 |          52.75 s |          0 |
| `regtest-gauxc-gapw-base` | 2 MPI x 2 OpenMP, CPU          |            3/3 |          16.90 s |          0 |
| `regtest-gauxc-cuda`      | 2 MPI x 2 OpenMP, two A40 GPUs |          22/22 |          50.96 s |          0 |

These sets cover direct valence, `PAW_ONE_CENTER`, `CP2K_DEFAULT`, `METHOD GAPW_XC`, MODEL NONE
isolation, the `PAW_ONE_CENTER_SPLIT` diagnostic, targeted force/stress diagnostics, inversion
reduction, the full internal k-point symmetry reduction, and the SPGLIB reduction. The three ECP
inversion, internal-symmetry, and SPGLIB inputs agree at `-13.36515054 Ha` in the final run. The
moved forward reconstruction also preserves the existing all-electron H2O XRD charge reference of
`-9.9999991403` electrons, and the existing mixed all-electron/GTH GAPW XRD input completes
successfully.

The original CPU k-point runs exposed a numerical-runtime defect rather than a SKALA or k-point
error. The packaged LibTorch uses oneMKL's grouped CBLAS ABI, while CP2K and ScaLAPACK use OpenBLAS;
LibTorch's calls to `cblas_sgemm_batch` and `cblas_dgemm_batch` were resolved to incompatible
OpenBLAS symbols and corrupted subsequent complex linear algebra. For LP64 OpenBLAS, CP2K now
expands these grouped operations into calls to the standard CBLAS GEMM interface. The all-electron
k-point test and all eight ECP tests then pass without `LD_PRELOAD`. Preloading an external oneMKL
remains invalid because it can interpose the complex BLAS used by ScaLAPACK.

The atom-composite CUDA path reuses the device storage of its coordinate, weight, density, gradient,
and kinetic-energy-density tensors when their shapes are unchanged. A paired H2/GTH `ENERGY_FORCE`
calculation from the same atomic density gives total energies of `-0.967804783962676 Ha` on the
uncached CPU path and `-0.967804781763126 Ha` on the cached CUDA path. Their difference is
`2.20e-9 Ha`; the largest difference among all six Cartesian force components is approximately
`3e-9 Ha/bohr`. Thus resetting the reusable leaves clears stale gradients and copies the new
primitive fields without changing the forward or backward result at chemically relevant precision.
The dedicated CUDA regression test exercises this molecular `PAW_ONE_CENTER` energy-and-force path.

The distributed-memory implementation evaluates one or more complete atom blocks on every active
rank. With `NATIVE_GRID_USE_CUDA T` and automatic device selection, MPI-local ranks select distinct
visible GPUs, so the same exact atom decomposition provides multi-GPU execution without gathering
model tensors to a root. A future refinement can accumulate the smooth PW interpolation adjoint
directly on the distributed PW layout instead of materializing that remaining full-grid object on
every rank.

The final two-A40 H2 check prints `0:0 1:1`, proving that both ranks execute on distinct devices.
Its one- and two-GPU energies are `-0.967804781763126` and `-0.967804783080038 Ha`; the difference
is `1.32e-9 Ha`, and the largest printed force-component difference is approximately `6e-9 Ha/bohr`.
The wall times are `6.79` and `6.74 s`. Peak memory is `2251 MiB` on the active device for one GPU
and `2339 MiB` per device for two GPUs, consistent with one model replica and one local atom block
per rank. All eight CUDA-regression outputs print `0:0 1:1`; across the complete 22-matcher set the
maximum sampled allocations are `11429` and `10553 MiB`, and both GPUs reach 100 percent sampled
utilization.

The original molecular atom-composite driver first evaluated SKALA on the complete regular grid and
then replaced that energy and potential with the atom-composite result. The first evaluation was
therefore redundant and made a 32-atom H test exceed one A40 after approximately `42.65 GiB` was
already allocated and another `4.97 GiB` was requested. The corrected driver prepares zeroed PW
potential storage and defers the only SKALA evaluation to the atom-composite path. The same test now
uses at most `4234 MiB` and takes `24.31 s` on one A40. With two MPI ranks and two A40s it uses
`2908` and `2032 MiB`, takes `16.09 s`, and obtains a `1.51`-fold speedup. The one- and two-GPU
energies differ by `6.03e-9 Ha`, and the largest difference over all 96 printed force components is
`4.0e-9 Ha/bohr`. Both outputs contain no regular-grid `skala_gpw_eval` timer.

All regular CPU validation uses the official Skala 1.1 checkpoint whose SHA-256 checksum is
`0c8432ac3f03c8f1276372df9aca5b7ee7f8939d47a8789eb158976e89aa0606`; the corresponding CUDA
checkpoint checksum is `f77be6002d873c0a2384b6df7850d32bbec519036344ff5fdde9730c6f9a4326`. The CPU
checksum matches the CP2K toolchain. The sparse atomic-weight taper intentionally changes the old
coarse 80-Ry k-point references. It simultaneously reduces the total translational-force residual
from the previous `10^-3 Ha/bohr` range to between `5.3e-8` and `1.5e-6 Ha/bohr`, which confirms
that the new references represent the continuous partition rather than a model-checkpoint change.

## Current acceptance status

- The 15 total gradient component-pair terms are generated by summing primitive fields before
  nonlinear feature construction. Exactly 14 of these pairs are new relative to the smooth
  total-gradient baseline.
- H2, H2O, and Ar atom-centered composite energies have been compared with fixed-density molecular
  GauXC all-electron references. H2 and Ar converge with grid refinement; H2O retains a documented
  finite GAPW reconstruction error.
- The composite representation preserves particle number algebraically. Its finite atom-grid
  integral is regression-tested explicitly; no density-dependent renormalization is introduced.
- Primitive-field directional derivatives and the complete PW plus one-center VXC adjoint agree with
  energy finite differences and matrix contractions.
- Molecular atom-composite forces pass finite differences for H2, H2O, Ar, and compressed H2, with
  the heavy-atom absolute error decreasing under grid refinement. The H2 molecular virial also
  passes an independent affine coordinate-scaling finite difference.
- Molecular `PAW_ONE_CENTER` with GTH-NLCC passes a bond-force finite difference, and its residual
  translation force decreases under grid refinement. Periodic native-grid PAW-like GAPW with NLCC
  passes targeted force and diagonal-stress finite differences. In both routes NLCC augments density
  and density gradient before feature construction while tau remains valence-only.
- The periodic regular-grid analytical virial passes full skew-cell finite-difference checks for
  both GPW and the all-electron GAPW common-grid composite after making sparse smooth-partition rows
  vanish continuously in SKALA's internal atom-grid integral. Periodic pseudopotential GAPW forces
  and stress also pass targeted GTH and ECP finite differences, with the ECP force residual
  decreasing by more than one order of magnitude under grid refinement. A periodic lattice-image
  construction for the separate atom-centered composite remains open.
- The periodic regular-grid GTH/ECP paths and their k-point symmetry reductions pass the production
  regression sets with two MPI ranks after removing an implicit `rho`/`rho_xc` mismatch and
  duplicate force/virial reductions. The molecular atom-centered composite gives rank-invariant
  energy and forces with one, two, and four MPI ranks. Complete atom blocks are evaluated on their
  owning CPU rank or GPU and device tensors are reused when their shapes are unchanged. Distributed
  PW-adjoint storage remains open and must preserve the validated energy, VXC, force, and virial
  behavior.
