# Nudged Elastic Band

The nudged elastic band (NEB) method finds a minimum-energy connection between two known structures.
In CP2K, select [RUN_TYPE BAND](#CP2K_INPUT.GLOBAL.RUN_TYPE) and configure
[MOTION/BAND](#CP2K_INPUT.MOTION.BAND). Intermediate replicas are optimized while springs control
their spacing. [BAND_TYPE](#CP2K_INPUT.MOTION.BAND.BAND_TYPE) selects, among other methods,
improved-tangent NEB (`IT-NEB`) or climbing-image NEB (`CI-NEB`). In the latter, the highest-energy
interior image climbs towards the saddle after [NSTEPS_IT](#CP2K_INPUT.MOTION.BAND.CI_NEB.NSTEPS_IT)
initial improved-tangent steps.

For fixed-cell Cartesian bands, [INTERPOLATION_METHOD](#CP2K_INPUT.MOTION.BAND.INTERPOLATION_METHOD)
selects how missing images are initialized. The default, `LI_BISECT`, bisects the longest interval.
Use `LI_IDPP` to relax the interpolated images on the image-dependent pair potential, or `S_IDPP` to
grow the band inward from the supplied images. Neither IDPP method requires physical forces during
initialization. For fixed-cell examples, see the tutorials below.

## Solid-state NEB

[SOLID_STATE](#CP2K_INPUT.MOTION.BAND.SOLID_STATE) extends the band to coupled atomic and cell
motion. Each replica has its own cell. At hydrostatic external pressure $p$, energies used for the
tangent and climbing-image selection are enthalpies $H=E+pV$, and analytical atomic forces and
stress provide their derivatives. This follows the joint atom/cell approach of
[Sheppard et al.](#Sheppard2012), using the fixed reference coordinates described below.

Use this mode for transformations between fully periodic crystals with the same composition and a
known correspondence between their atoms and lattice vectors. It does not identify that
correspondence automatically. Optimize the two endpoints separately at the intended pressure before
constructing the band. Both endpoint positions and cells remain fixed during the band optimization.

### Preparing the input

In addition to a force-evaluation method and `RUN_TYPE BAND`, set the following:

- `FORCE_EVAL/STRESS_TENSOR ANALYTICAL` for the full stress tensor. The current implementation
  supports Quickstep and FIST.
- `SOLID_STATE TRUE`, `POT_TYPE FULL`, and at least three replicas.
- `BAND_TYPE IT-NEB` or `CI-NEB`, with `OPTIMIZE_BAND/OPT_TYPE DIIS` and
  `OPTIMIZE_END_POINTS FALSE`.
- [EXTERNAL_PRESSURE](#CP2K_INPUT.MOTION.BAND.EXTERNAL_PRESSURE) explicitly, in bar by default. Set
  zero for the potential-energy surface without an external pressure contribution.
- `PERIODIC XYZ` for `FORCE_EVAL/SUBSYS/CELL`. Do not request cell symmetry or `CANONICALIZE TRUE`.
  For explicitly oriented cells, `CANONICALIZE FALSE` preserves the supplied orientation.

Every input [REPLICA](#CP2K_INPUT.MOTION.BAND.REPLICA) needs
[CELL_MATRIX](#CP2K_INPUT.MOTION.BAND.REPLICA.CELL_MATRIX). Its nine entries are the Cartesian
components of vectors A, B, and C, in that order, in angstrom by default. For example,
`CELL_MATRIX 10 0 0 0.8 9.7 0 -0.3 0.4 10.2` specifies $\mathbf A=(10,0,0)$,
$\mathbf B=(0.8,9.7,0)$, and $\mathbf C=(-0.3,0.4,10.2)$. Replica `COORD` entries contain only
Cartesian coordinates, without element or kind labels, in the same atom order as
`FORCE_EVAL/SUBSYS/COORD`. Alternatively, use `COORD_FILE_NAME` with an XYZ file. The cell must
still be supplied through `CELL_MATRIX`.

Provide right-handed endpoint cells in a common orientation and with the same lattice
correspondence. For example, both can have A along X and B in the XY plane. A common rigid rotation
of all cells and positions is also allowed. Do not rotate the endpoints independently. Use
consistent, unwrapped atomic coordinates across the band, rather than wrapping the atoms of each
replica independently into its own cell.

At least two input replicas are required. If
[NUMBER_OF_REPLICA](#CP2K_INPUT.MOTION.BAND.NUMBER_OF_REPLICA) is larger than the number supplied,
CP2K inserts replicas by bisection in the joint atomic and cell coordinates. This interpolates
logarithmic cell stretches, not Cartesian positions and raw cell matrices independently. Inspect the
resulting structures for close contacts and inappropriate atom mapping.

### Production input for a variable-cell band

For a pair of endpoints relaxed at zero external pressure, the following `MOTION` section starts
with seven replicas and turns on the climbing image after ten improved-tangent iterations. Supply
`initial.xyz` and `final.xyz` in matching atom order. Replace the example cell matrices with the
respective relaxed cells. Use the same force-evaluation settings as in the endpoint calculations,
including analytical stress.

```none
&MOTION
  &BAND
    SOLID_STATE TRUE
    BAND_TYPE CI-NEB
    NUMBER_OF_REPLICA 7
    EXTERNAL_PRESSURE [bar] 0
    CELL_WEIGHT 1.0
    &CI_NEB
      NSTEPS_IT 10
    &END CI_NEB
    &OPTIMIZE_BAND
      OPT_TYPE DIIS
      OPTIMIZE_END_POINTS FALSE
      &DIIS
        MAX_STEPS 200
      &END DIIS
    &END OPTIMIZE_BAND
    &CONVERGENCE_CONTROL
      MAX_FORCE 4.5E-4
      RMS_FORCE 3.0E-4
      MAX_DR 2.0E-4
      RMS_DR 1.0E-4
    &END CONVERGENCE_CONTROL
    &REPLICA
      COORD_FILE_NAME initial.xyz
      CELL_MATRIX 10 0 0 0 10 0 0 0 10
    &END REPLICA
    &REPLICA
      COORD_FILE_NAME final.xyz
      CELL_MATRIX 10.5 0 0 0.8 9.7 0 -0.3 0.4 10.2
    &END REPLICA
  &END BAND
&END MOTION
```

Check the reported convergence of the joint atomic and cell coordinates. Increase the replica count
and tighten the force thresholds until the activation enthalpy is stable at the required accuracy.

For Quickstep, keep the basis, functional, occupations, k-point sampling and other electronic
settings consistent along the band. A suitable
[CELL_REF](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.CELL_REF) can keep the plane-wave grid dimensions
fixed while cells change. Converge energies, forces and stresses with respect to the grid and SCF
tolerance. Finite-basis and grid errors in stress can otherwise change the band.

### Coordinates, metric and convergence

The first replica defines a fixed cell frame. With cell vectors as columns, write its cell matrix as
$h_0=Q U_0$, where $Q$ is a fixed rotation and $U_0$ is upper triangular. Every replica is
represented by

$$
h=Q D U_0,\qquad \mathbf r_i=Q D\mathbf u_i ,
$$

where $D$ is upper triangular with positive diagonal entries. The atomic coordinates
$\mathbf u_i=U_0 h^{-1}\mathbf r_i$ are proportional to fractional coordinates. The six cell
coordinates are

$$
\mathbf q=L_c\left(\ln D_{11},\ln D_{22},\ln D_{33},D_{12},D_{13},D_{23}\right),\qquad
L_c=w\sqrt{N}\left(\frac{V_0}{N}\right)^{1/3}.
$$

Here $N$ is the number of atoms, $V_0=\det h_0$, and $w$ is
[CELL_WEIGHT](#CP2K_INPUT.MOTION.BAND.CELL_WEIGHT), which defaults to one. The orientation,
reference cell and length scale remain fixed throughout the optimization. Logarithmic diagonal
coordinates prevent changes of handedness without introducing three redundant cell rotations.

The band uses Euclidean distances in the concatenated vector
$(\mathbf u_1,\ldots,\mathbf u_N,\mathbf q)$. The same metric enters tangents, springs, force
projection and climbing-image updates. `CELL_WEIGHT` changes the relative cost of cell and atomic
motion, and can therefore change the resulting connection. It is not just an optimizer step size. It
should be reported with the calculation and kept unchanged on restart.

Forces are the negative derivatives of $E+pV$ with respect to these coordinates. The configurational
virial already includes affine atomic motion at fixed fractional coordinates, so no additional
atomic force-position contribution is added to the cell derivative.

The displacement and force tolerances in
[CONVERGENCE_CONTROL](#CP2K_INPUT.MOTION.BAND.CONVERGENCE_CONTROL) apply to this joint vector. Their
values do not directly specify a Cartesian displacement or a pressure tolerance. Check the residual
band forces as well as whether the optimization actually converged.

### Output and restarts

The replica information reports the potential energy, enthalpy, cell volume and six cell-force
components. Solid-state band summaries and the band energy file use enthalpies, including the
external pressure contribution. A band total also includes the artificial spring energy and is not
the activation enthalpy. Obtain forward and reverse barriers from the saddle enthalpy relative to
the respective endpoint.

[FINAL_BAND](#CP2K_INPUT.MOTION.BAND.FINAL_BAND) writes extended XYZ with a separate `Lattice` and
`Enthalpy` for every replica. The motion restart stores every replica's coordinates and
`CELL_MATRIX`. Restart from that complete band state, not just the last cell printed by a
force-evaluation process. Keep the first endpoint, pressure and `CELL_WEIGHT` unchanged. The DIIS
history is not stored, so the restarted optimization need not reproduce the identical sequence of
steps.

### Limitations

The initial implementation requires full periodicity, analytical stress, DIIS, IT-NEB or CI-NEB,
fixed endpoints, and an unconstrained atom-only system. It does not support shell particles,
colvar-projected bands, numerical or diagonal-only stress, nonhydrostatic external stress, or
automatic matching of atoms and inequivalent supercells. Neither `LI_IDPP` nor `S_IDPP` can be
combined with `SOLID_STATE`.

Atomic `ALIGN_FRAMES` and `ROTATE_FRAMES` are disabled in solid-state mode. The common cell frame
determines the orientation, and rotating atomic coordinates alone would change the crystal.

## Tutorials

- <https://www.cp2k.org/exercises:common:neb>
- <https://www.cp2k.org/exercises:2018_uzh_cmest:path_optimization_neb>
- <https://www.cp2k.org/exercises:2017_ethz_mmm:nudged_elastic_band>
- <https://www.cp2k.org/exercises:2015_cecam_tutorial:neb>

```{youtube} t8jxq4e0Qtw
---
url_parameters: ?start=2
align: center
privacy_mode:
---
```
