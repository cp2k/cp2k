# Density of States

The density of states (DOS) counts how many Kohn-Sham states occur at a given energy. In practice it
is used together with band structures, occupations, and orbital projections to identify the
electronic states responsible for a material or molecular property.

In electronic-structure analysis, DOS and projected DOS (PDOS) are most often used to answer
questions such as:

- Is the system metallic, semiconducting, or insulating? A finite DOS at the Fermi level usually
  indicates metallic behavior, while a gap around the band edges indicates an insulator or
  semiconductor.
- Which atoms and orbitals form the valence-band maximum, conduction-band minimum, frontier
  molecular orbitals, or states near the Fermi level?
- Are defect, dopant, adsorbate, or surface states inside the gap, and are they localized on
  particular atoms?
- Which angular-momentum channels dominate a feature, for example transition-metal `d` states,
  oxygen `p` states, or adsorbate frontier orbitals?
- How strong is orbital hybridization, such as metal--oxygen `d`--`p` mixing, substrate--adsorbate
  coupling, or ligand-field splitting?
- How do strain, oxidation state, magnetic ordering, charging, adsorption, or structural relaxation
  shift band edges and redistribute orbital character?

DOS gives the total number of available states. PDOS decomposes the same spectrum into selected
atomic kinds, angular-momentum channels, individual components, atom lists, or real-space regions.
For this reason, total DOS is usually used to locate the energy windows of interest, while PDOS is
used to assign the chemical character of the peaks or band edges.

```{important}
Most of the DOS/PDOS functionality described below, including the unified
`&DFT%PRINT%DOS` interface and the `&CURVE` subsection for broadened, directly plottable
DOS/PDOS curves, is only available in CP2K 2026.2 and later versions. The ordinary `.dos` and
`.pdos` outputs keep the traditional non-broadened data, while `&CURVE` requests additional
broadened curve output.
```

## Basic CP2K input

For DFT calculations, DOS output is controlled by
[&DFT%PRINT%DOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS). Requesting this print section writes the
total DOS. Set [&PDOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.PDOS) under the same section to
additionally write projected DOS files.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &DOS
        DELTA_E [eV] 0.01
        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO AUTO
          &BROADEN
            WIDTH [eV] 0.1
          &END BROADEN
        &END CURVE
        &PDOS
          COMPONENTS F
        &END PDOS
      &END DOS
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

Useful input keywords are documented in the input reference:

- [DELTA_E](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.DELTA_E): choose the energy-grid spacing.
- [NLUMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.NLUMO): requested number of unoccupied states for
  DOS/PDOS.
- [&CURVE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.CURVE): request additional broadened, directly
  plottable DOS/PDOS curve output.
- [ENERGY_UNIT](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.CURVE.ENERGY_UNIT): print curve energies in
  Hartree or eV.
- [ENERGY_ZERO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.CURVE.ENERGY_ZERO): choose the reference energy
  for curve output.
- [&BROADEN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.CURVE.BROADEN): choose the line shape and
  broadening width for curve output.
- [&PDOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.PDOS): request projected DOS on each element.
- [&LDOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.LDOS) and
  [&R_LDOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.R_LDOS): request atom-list or real-space local
  projected DOS.

For PDOS, there are some additional controls:

- [COMPONENTS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.PDOS.COMPONENTS): split angular-momentum
  channels into individual components. This keyword also applies to `&LDOS` and `&R_LDOS`.

```{note}
For diagonalization-based SCF calculations, the unoccupied states used for DOS/PDOS are those
available from the Kohn-Sham diagonalization, controlled by
[ADDED_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.ADDED_MOS). In this case, `NLUMO` acts mainly as a lower
bound: CP2K may increase the number of additional MOs if needed, but not truncate unoccupied MOs
that are already available. In OT calculations, unoccupied states for DOS/PDOS are computed in a
separate post-SCF step.
```

### Choosing the energy reference

Use [ENERGY_ZERO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.CURVE.ENERGY_ZERO) to make different curve
plots comparable. This setting applies to broadened curve output, while ordinary non-curve outputs
are not shifted by this setting.

- `FERMI` is normally the appropriate zero for metals and smeared calculations.
- `HOCO` places the highest occupied crystal orbital at zero. It is often more convenient for
  molecules, semiconductors, and insulators, especially when comparing DOS with band plots where the
  valence-band maximum is set to 0 eV.
- `ABSOLUTE` prints the absolute Kohn-Sham eigenvalues.
- `AUTO` selects `FERMI` when smearing or fractional occupations are present, and `HOCO` otherwise.

The curve output header reports both reference energies and the selected zero, for example:

```text
# E(Fermi) = 0.223578 a.u. = 6.08391 eV
# E(HOCO)  = 0.190000 a.u. = 5.17018 eV
# Energy zero: AUTO -> HOCO
```

When comparing several structures or charge states, choose the same energy-zero convention for all
curve plots. In spin-polarized calculations, the two spin channels use the same reference energy so
that alpha and beta curves remain on a common energy axis. If an absolute alignment is needed across
different cells or surfaces, additional electrostatic-potential alignment may be required; changing
the DOS zero alone is not a substitute for such an alignment.

## Broadening, k-points, and gaps

Electronic-structure features observed in DOS/PDOS, such as peak positions, band gaps, and peak
widths, can be affected by the calculation setup. Appropriate choices of the functional, basis set,
k-point sampling, smearing, unoccupied states, and broadening parameters are important for obtaining
a reliable interpretation of the electronic structure.

The broadened curve output replaces each discrete eigenvalue by a normalized line shape. The
broadening [WIDTH](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.CURVE.BROADEN.WIDTH) is the full width at
half maximum (FWHM). The default value of 0.1 eV is a reasonable choice for most applications.
Smaller values can be useful for resolving narrow features, small gaps, spin splittings, or defect
states, whereas larger values mainly produce smoother curves and may hide such details. Curve output
requires `DELTA_E > 0`; if `DELTA_E <= 0`, CP2K skips the broadened curve output and prints a
warning.

DOS and PDOS usually require a denser k-point mesh than calculations aimed only at total energies,
forces, or geometry optimization. A too sparse k-point mesh gives a DOS dominated by the discrete
sampling of the Brillouin zone, especially near the band edges or the Fermi level. Broadening helps
produce a smooth curve, but it should not be used as a substitute for insufficient k-point sampling.

For molecular or very large supercell calculations with only the Gamma point, the broadened curve is
a representation of discrete levels. It can still be useful for assigning orbital character, but
peak widths are chosen mainly for visualization rather than for Brillouin-zone integration.

## Interpreting PDOS

With `&PDOS`, CP2K writes projected DOS, usually one file per atomic kind and spin channel. The
ordinary `.pdos` files contain state-resolved projected weights. If `&CURVE` is also requested, CP2K
additionally writes broadened PDOS curve files with a projected `total` column followed by
angular-momentum or component-resolved columns:

```text
# Energy[eV]  total  s  p  d
```

With [`COMPONENTS T`](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.PDOS.COMPONENTS), the channels are split
into individual components:

```text
# Energy[eV]  total  s  py  pz  px  ...
```

The `total` column in one PDOS curve file is the projected total for that atomic kind, not the total
DOS of the full system. Summing the projected totals over all kinds should reproduce the total DOS
up to numerical and output precision, provided that the same states and broadening are used.

```{warning}
Component-resolved k-point PDOS with full k-point symmetry reduction should be interpreted with care,
because symmetry-reduced k-points can rotate the orbital components. It is suggested to disable
k-point symmetry when analyzing individual components.
```

Practical PDOS analysis usually proceeds by selecting an energy window and assigning the dominant
contributions:

- band-edge states: compare the PDOS around the valence-band maximum and conduction-band minimum;
- defect or dopant levels: check whether in-gap peaks are localized on the defect atoms or nearby
  host atoms;
- surface and adsorption states: compare slab atoms, adsorbate atoms, and bulk-like atoms;
- transition-metal compounds: inspect the relative positions and spin splitting of metal `d` and
  ligand `p` channels;
- molecular frontier levels: identify which fragments contribute to HOMO/HOCO and LUMO/LUCO
  features.

For quantitative orbital populations, integrated PDOS values should be used with care because they
depend on the projection scheme, basis set, and chosen energy window.

State-resolved k-point PDOS output is currently not available. K-point PDOS is accumulated with
k-point weights for the broadened curve output.

### Projection used for PDOS

CP2K's standard PDOS is a Löwdin-type projection in the atomic-orbital basis. For non-k-point
calculations, the molecular-orbital coefficients `C` are first transformed with the square root of
the overlap matrix `S`,

$$
\widetilde C = S^{\frac{1}{2}} C
$$

and the projected weight of a state is accumulated from $|\widetilde C_{\mu i}|^2$ over atomic
orbitals belonging to a given atomic kind and angular-momentum channel. This gives a positive,
orthogonalized AO decomposition of the DOS. It is not a Mulliken population analysis.

For k-point calculations, CP2K performs the same Löwdin-type projection separately at each k point,
using the corresponding overlap matrix and Bloch MO coefficients. The projected contributions are
then accumulated with the k-point weights to form the broadened PDOS.

## Plotting DOS/PDOS

### Using gnuplot

The broadened DOS and PDOS curve files can be plotted directly with [gnuplot](http://gnuplot.info/),
since comment lines starting with `#` are ignored. For example, if the first column is the energy
and the second column is the total DOS, a simple plot is

```gnuplot
plot "project-curve-1.dos" using 1:2 with lines title "DOS"
```

The plotted energy range can be restricted directly in the `plot` command. For example, to show only
states within 5 eV of the chosen energy zero,

```gnuplot
plot [-5:5] "project-curve-1.dos" using 1:2 with lines title "DOS"
```

Both the energy and DOS ranges can be limited in the same way:

```gnuplot
plot [-5:5][0:*] "project-curve-1.dos" using 1:2 with lines title "DOS"
```

For PDOS curve files, different columns can be plotted together. For example, if the second column
is the kind-resolved total PDOS and the following columns are angular-momentum components,

```gnuplot
plot [-5:5] "project-k1-curve-1.pdos" using 1:2 with lines title "total", \
            "" using 1:3 with lines title "s", \
            "" using 1:4 with lines title "p"
```

When comparing DOS/PDOS curves with band structures, use the same energy zero, for example `FERMI`
for metals or `HOCO` for many insulating and molecular systems.

## Local projected DOS

[&LDOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.LDOS) projects onto explicitly listed atoms.
[&R_LDOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS.R_LDOS) projects onto a real-space volume defined
relative to listed atoms. These outputs are useful for checking whether a DOS feature is localized
in a defect region, adsorbate, surface layer, or selected fragment.They are currently not
implemented for k-points.

```{caution}
The definition of `LDOS` in CP2K is different from the real-space local density of states
$\rho(\mathbf{r}, E)$ often used in, for example, STM-related analysis. In this DOS/PDOS print
section, `LDOS` denotes a projected DOS for user-defined lists of atoms, i.e. a local PDOS over
selected atoms. `R_LDOS` further restricts the analysis to a user-defined spatial region around
selected atoms. These outputs are integrated projected quantities, not three-dimensional
real-space LDOS maps.
```

## XAS and XAS_TDP projected DOS

XAS and XAS_TDP have their own dedicated [&PDOS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS.PRINT.PDOS) print
sections. They are intentionally separate from the DFT [&DOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.DOS)
section because their output is used for spectral-feature analysis rather than ordinary ground-state
density-of-states analysis.

For XAS calculations based on transition-potential or core-hole electronic structures, the projected
spectrum is computed from the XAS electronic structure. It can be used to inspect which atoms and
angular-momentum channels contribute to states relevant for the near-edge spectrum. The
interpretation is therefore tied to the core-excited or transition-potential state, not to the
ground-state DOS.

For XAS_TDP, the distinction is stronger. The `.pdos` files reuse the standard PDOS projection and
file format, but the states being projected are linear-response orbitals and the energy column
represents excitation energies. Occupation numbers and references to molecular orbitals in the
generic PDOS format should not be interpreted as ordinary ground-state Kohn-Sham occupations.

A typical XAS_TDP use is to assign intense spectral features to atomic or orbital character, for
example to distinguish metal-centered, ligand-centered, or mixed transitions. In this context,
`DELTA_E` and `&CURVE` remain useful plotting controls, while ordinary DFT-DOS concepts such as
Fermi-level alignment, HOCO alignment, or total DOS are usually not the relevant analysis target.
See [&XAS%PRINT%PDOS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS.PRINT.PDOS) and
[&XAS_TDP%PRINT%PDOS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.PRINT.PDOS) for the corresponding input
sections.
