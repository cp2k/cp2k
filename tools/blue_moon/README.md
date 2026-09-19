# Blue-moon postprocessing for one collective constraint

This standalone NumPy tool addresses the common single-coordinate cases in
[CP2K issue #5863](https://github.com/cp2k/cp2k/issues/5863), without changing the MD engine. It
computes the scalar mass metric and its derivative from the saved geometry, then combines them with
CP2K's SHAKE multiplier:

```text
g = grad(xi), H = Hessian(xi), v = M^-1 g
Z = g . v
G = v . H . v / Z^2
w = Z^(-1/2)
dA/dxi = sum[w * (-lambda_SHAKE + kB*T*G)] / sum[w]
```

The derivatives are obtained by second-order forward automatic differentiation of the coordinate
definition; no finite-difference step size is needed. All occurrences of a shared atom use the same
Cartesian independent variables. Plane atoms and bond-center atoms move in the differentiation too.
The implementation is intentionally separate from CP2K and can be adapted or incorporated into a
postprocessing package.

## Supported scope and prerequisites

- Exactly **one** fixed collective constraint in an equilibrated, fixed-cell trajectory. Additional
  fixed atoms, rigid molecules, other SHAKE constraints, PIMD, RESPA, moving targets and changing
  cells are **not supported**. The log parser rejects multiple multipliers, but cannot detect fixed
  atoms or other omitted simulation settings: check the original CP2K input. This tool does not
  parse or validate that input.
- Standard velocity-Verlet output with one SHAKE/RATTLE pair per recorded step. RATTLE is validated
  but never counted as another configurational-force sample.
- A CP2K XYZ position trajectory in **Angstrom**, retaining the `i = ...,` step labels. The atomic
  order, geometry representation, cell and image choices must match the run. Use the original
  trajectory, not a wrapped/reordered visualization export. In particular, reconstructing molecular
  images can change centroid-based CVs.
- The **actual masses used by CP2K**, including isotope or KIND/MASS overrides, for every
  participating real atom. All must be finite and positive; element labels are not used to guess
  masses. Common conversion from amu to electron masses cancels in G and in the normalized weighted
  average, so Z is reported with masses in amu.
- The equilibrium temperature in kelvin, and the exact fixed CV target and tolerance in CP2K
  internal coordinate units (bohr for distances, radians for angles). Mixed-coordinate coefficients
  must reproduce the units and normalization of the original `COMBINE_COLVAR`. Changing the
  normalization without transforming lambda changes the result.

## Run

Requires Python 3.9+ and NumPy 1.21+ (including `numpy.typing`). From the repository root:

```shell
python3 tools/blue_moon/blue_moon.py \
  tools/blue_moon/examples/distance.json \
  tools/blue_moon/examples/trajectory.xyz \
  tools/blue_moon/examples/constraint.LagrangeMultLog \
  --first-step 1
```

The included five-step, two-argon-atom CP2K 2026.1 output is only an I/O smoke test, **not an
equilibrated free-energy calculation**. Its generating input is `examples/distance.inp`.
`examples/coordinates.inp` independently checks all six coordinate forms against CP2K's
`METADYN/COLVAR` output, without depositing bias hills. Run those inputs in a scratch directory. A
suitable production configuration looks like:

```json
{
  "cv": {
    "type": "linear_combination",
    "terms": [
      {"coefficient": 1, "cv": {"type": "distance", "atoms": [1, 2]}},
      {"coefficient": -1, "cv": {"type": "distance", "atoms": [3, 2]}}
    ]
  },
  "cell_angstrom": [[20, 0, 0], [0, 20, 0], [0, 0, 20]],
  "masses_amu": {"1": 15.9994, "2": 1.00794, "3": 15.9994},
  "temperature_kelvin": 300,
  "target_au": -1.0,
  "target_tolerance_au": 0.00001
}
```

Substitute the values from your simulation. Atom indices are one-based and masses are keyed by those
indices. Cell vectors A, B, C are **rows**, in Angstrom. A fixed nonsingular cell must be provided
even for isolated systems, since CP2K's DISTANCE, ANGLE and TORSION implementations apply
minimum-image wrapping using the cell. The tolerance must accommodate constraint and trajectory
output precision, while remaining tight enough to detect a mismatched CV.

The multiplier file has **no step labels**. `--first-step` is therefore mandatory: it specifies the
MD step of the first SHAKE record, not the first XYZ frame. For a fresh standard run without
`CONSTRAINT_INIT`, this is normally 1; the initial XYZ frame at step 0 is ignored. Do not infer this
offset from matching array lengths. Inspect initialization/restart behavior and split appended runs
before analysis. Remove any initialization constraint records explicitly.

Print XYZ and multipliers at the same cadence (preferably every MD step). Set `--stride` to that
cadence if it differs from 1. Every XYZ frame at or after `--first-step` must have a matching pair;
gaps, count mismatches and repeated/decreasing XYZ steps are rejected. Initial XYZ frames before
`--first-step` are ignored, but no SHAKE records are silently dropped. `--discard N` discards N
**paired samples** as equilibration, after alignment and validation.

CSV on standard output contains step, CV, lambda, Z, G, weight and instantaneous corrected force.
The final JSON summary is on standard error and is printed **only on success**. On failure, discard
any partial CSV and inspect the nonzero exit status. Free-energy gradients are in hartree per
internal CV unit; for a distance, divide by `0.52917720859` to convert from hartree/bohr to
hartree/Angstrom. This program does not integrate windows or claim statistical convergence. Check
time-step/constraint-tolerance convergence and use block analysis or block bootstrap of the weighted
numerator **and** denominator for correlated uncertainty estimates.

## Coordinate definitions

| `type`               | `atoms` ordering and CP2K correspondence                                                                                  |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `distance`           | `[i,j]`: ordinary, unsigned `DISTANCE` with `AXIS XYZ`; also usable as a term in `DISTANCE_FUNCTION`                      |
| `angle`              | `[i,j,k]`: `ANGLE`, with j at the vertex, in radians                                                                      |
| `torsion`            | `[a,b,c,d]`: `TORSION`, with CP2K's sign; also requires `"reference": angle_in_radians` to choose its continuous branch   |
| `point_plane`        | `[i,j,k,l]`: `DISTANCE_POINT_PLANE`, plane atoms i,j,k followed by point l; signed normal `(ri-rj) x (rk-rj)`             |
| `point_bond_center`  | `[p,i,j]`: `DISTANCE` from atom p to a `POINT TYPE GEO_CENTER` of i,j; this is an arithmetic, not mass-weighted, midpoint |
| `linear_combination` | `terms` containing `coefficient` and another `cv`; any number of terms, including nested combinations                     |

This covers the eight common coordinate forms listed in the issue, **not** arbitrary CP2K
`COMBINE_COLVAR` expressions. No Python expression evaluation or automatic CP2K input translation is
performed. Unsupported types, options and keys are rejected.

Each primitive defaults to `"pbc": true`, using CP2K's fractional-coordinate minimum-image
convention, including triclinic cells. `"pbc": false` is available for matching
`DISTANCE_FUNCTION/PBC FALSE` terms or `DISTANCE_POINT_PLANE/PBC FALSE`; do not use it to disable
wrapping that the original CP2K coordinate applies. Half-cell branch boundaries, zero distances and
degenerate angle/plane/torsion geometries are rejected. A torsion `reference` must select the same
branch as CP2K throughout the trajectory; derivatives do not fix a wrong branch inside a
combination. The fixed-target check provides an additional consistency test.

For a point-plane CV the centroid displacement and the two plane vectors are wrapped separately, as
in CP2K. For the bond center the arithmetic mean is formed **before** wrapping its distance to the
point. Keep the bond constituents in the same molecular image as in the run. Other virtual-point
definitions, signed/projected distances and per-axis PBC are not supported. For torsions, the
minimum images of the 1-3 displacements must agree with the sum of adjacent bond images, as required
for CP2K's implemented gradient to match the differentiated coordinate; inconsistent images are
rejected.

## Validation and the distance-difference special case

```shell
python3 -m unittest discover -s tools/blue_moon -p 'test_*.py' -v
python3 -m mypy --strict tools/blue_moon/
```

Tests cover analytic distance and bond-center metrics, shared-atom distance differences,
finite-difference checks of every supported primitive's gradient and Hessian, mass/CV scaling,
triclinic images, torsion branches, the weighted estimator and strict input alignment. The
miscellaneous CI job runs both the unit tests and strict type checking for the complete directory,
including tests. Runtime validation of JSON input remains necessary; the numerical evaluator uses a
typed coordinate tree only after the input has been checked.

### Recorded CP2K distance-difference trajectory

`examples/distance_difference.inp` generates a 32-step NVE reference for three noninteracting
particles with explicit KIND masses of 12, 1 and 16 amu. It constrains `COMBINE_COLVAR R1-R2` to -1
Angstrom and also prints the individual distances, their angle and the equivalent
`DISTANCE_FUNCTION` without depositing metadynamics hills. The supplied initial velocities are
tangent to the constraint and have zero total momentum; `SHAKE_TOLERANCE` is tightened to `1e-12`.

The accompanying `.xyz`, `.LagrangeMultLog` and `.metadynLog` files were generated with **CP2K
2026.1, revision 5e54ba2**, one MPI process and one OpenMP thread. Reproduce them in an empty
scratch directory with:

```shell
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 cp2k.psmp -i /path/to/examples/distance_difference.inp -o reference.out
```

From the repository root, analyze the recorded data with:

```shell
python3 tools/blue_moon/blue_moon.py \
  tools/blue_moon/examples/distance_difference.json \
  tools/blue_moon/examples/distance_difference.xyz \
  tools/blue_moon/examples/distance_difference.LagrangeMultLog \
  --first-step 1 --discard 4
```

The regression checks every frame against the native CV output (allowing for its five printed
decimal places), independently evaluates Z and G using the distance-difference formulas below, and
tests the actual CLI's pairing, discard, units, CSV columns and final weighted reduction. Both Z and
G vary and G is nonzero. An initial-time inertial-force calculation using the specified velocities
also checks the SHAKE multiplier's sign and normalization within timestep and printed precision.
Regenerating with another CP2K version should be compared with these invariant checks; the recorded
data are not an assertion of bitwise reproducibility across versions and platforms.

For the initial-time check, reducing the timestep from 0.25 to 0.125 and 0.0625 fs changes the
printed first multiplier from `-3.092e-6` to `-3.091e-6` and `-3.091e-6` hartree/bohr, approaching
the continuous-time value `-3.0905654669e-6`. The nine-decimal-place multiplier output limits this
comparison. This checks the initial inertial force, not an equilibrated mean force.

For the recorded files, 28 retained samples give `2.21025541695497e-5` hartree/bohr versus the
uncorrected `3.2065e-6`. These numbers are **I/O/estimator regression values, not a physical
free-energy result**: this short unthermostatted, unconfined trajectory is not an equilibrium sample
at 300 K. The temperature and discarded prefix exercise the processing options only.

### Independent canonical-ensemble reference

`test_ensemble.py` separately checks the thermodynamic mean-force identity for a confined,
three-particle model with `U = k*(r12^2+r32^2)/2`, `xi = r12-r32`, and `k = kB*T/bohr^2`. In atomic
length units, removing translation/rotation constants gives

```text
P(xi) = integral_0^infinity r1^2 r2^2 exp(-(r1^2+r2^2)/2) ds
r1 = s + max(xi, 0), r2 = s + max(-xi, 0)
dA/dxi = -kB*T * d(ln P)/dxi
```

The reference integrates this one-dimensional expression by Gauss-Legendre quadrature and
finite-differences `ln P`. It contains neither Z nor G nor calls to the postprocessor. Quadrature
order and finite-difference step are checked separately.

For comparison, the constrained configurational measure is proportional to
`sqrt(Z)*r1^2*r2^2*exp(-U/(kB*T)) ds d(cos(theta))`. The test integrates the tangent-space Maxwell
velocities analytically. With `g = grad(xi)`, their covariance is
`C = kB*T*(M^-1 - (M^-1*g)*(M^-1*g)^T/Z)`. Twice differentiating the fixed constraint yields

```text
<-lambda | positions> = (g . M^-1 . grad(U) - trace(H*C))/Z
```

The gradients and radial Hessian blocks used here are constructed explicitly, independently of the
tool's automatic differentiation and metric functions. Those conditional multipliers and geometries
are passed through `analyze`; numerical integration of its weighted numerator and denominator is
then compared with `-kB*T*d(ln P)/dxi`.

The two windows `xi = -1` and `0.75` bohr give respectively `-1.1238227746` and `0.8527258234` for
`(dA/dxi)/(kB*T)` in inverse bohr. Both mass sets `[12,1,16]` and `[1,2,3]` recover the same
thermodynamic answer to within the test tolerance of `2e-7` inverse bohr, although their constrained
distributions differ. Removing either the G term or the Z reweighting fails this comparison.

This is deterministic equilibrium integration, **not a synthetic CP2K trajectory or a CP2K MD
convergence benchmark**. Together with the recorded-data test it checks complementary parts of the
workflow, but does not establish long-time sampling, thermostat or timestep convergence in CP2K.
Independent scientific review of the distance-difference specialization remains appropriate; these
tests do not substitute for that review.

### Equivalent logarithmic derivative

For one constraint the correction can also be expressed as `G = (1/2) D_xi ln Z`, provided the
derivative direction is specified:

```text
b = M^-1 grad(xi) / Z
D_xi = b . grad
D_xi xi = 1
G = (1/2) D_xi ln Z = b . grad(Z) / (2*Z)
```

This is the local mass-weighted normal derivative at each saved configuration, not a time derivative
along a constrained trajectory. Within a fixed window `xi` is constant, while `Z` can change with
other degrees of freedom (for example the angle between bonds at fixed distance difference). Ratios
of frame-to-frame changes `Delta ln Z / Delta xi`, or differences of window averages, therefore do
not in general give the required correction. Constraint-tolerance noise must not be used as the
denominator of that ratio.

The full Cartesian Hessian need not be stored to evaluate this identity: a Hessian-vector product or
a directional derivative of `Z` is sufficient. Such a derivative still requires differentiating the
coordinate gradients or evaluating them at displaced geometries; it is not determined by the values
of `Z` on the saved trajectory alone. This implementation retains the exact second-order
forward-differentiation formulation for the supported small coordinate definitions. An additional
regression compares its result with `(ln Z(r+h*b) - ln Z(r-h*b)) / (4*h)` at several step sizes for
all supported primitive types, combinations and a rescaled coordinate. It also verifies
`D_xi xi = 1`.

### Distance difference

For `xi = rij - rkj`, with `c = rho_ij . rho_kj`, direct differentiation gives:

```text
Z = 1/mi + 1/mk + 2*(1-c)/mj
G = (1-c*c) * (1/rij - 1/rkj) / (mj*mj*Z*Z)
```

G is not generally zero: equal distances or collinear bonds are special cases. For example,
positions `[(2,0,0),(0,0,0),(1,3,0)]` in bohr and masses `[12,1,16]` in amu give
`Z = 1.5133778012996573` and `G = 0.07221504489510591` per bohr. The tests compare the Hessian
contraction, this closed form and an independent directional finite difference of Z.

The general estimator is Eq. (6)-(8) of
[Komeiji (2007), Chem-Bio Informatics Journal 7, 12](https://doi.org/10.1273/cbij.7.12). The `G=0`
claim in that paper's Eq. (35) does not follow from Eq. (8) for a general three-atom distance
difference. The implementation uses the general expression, also consistent with the
single-constraint specialization of the
[VASP blue-moon expression](https://vasp.at/wiki/Blue_moon_ensemble), accounting for CP2K's opposite
lambda sign convention. See also the CP2K manual's
[constrained-dynamics chapter](https://manual.cp2k.org/trunk/methods/sampling/constrained_dynamics.html).
