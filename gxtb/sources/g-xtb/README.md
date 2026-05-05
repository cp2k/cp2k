# 🚧 g-xTB — Development Version

This is a preliminary version of g-xTB, a general-purpose semiempirical quantum mechanical method approximating ωB97M-V/def2-TZVPPD properties for all elements up to lawrencium (Z = 1–103).

## 📄 Preprint

See the [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/685434533ba0887c335fc974) preprint for the first g-xTB version. The present version includes substantial improvements and changes beyond the preprint.

## 📦 Installation

We provide statically linked binaries of the [xtb](https://github.com/grimme-lab/xtb) program ([modified version 6.7.1](https://github.com/thfroitzheim/xtb/tree/gxtb)), interfacing with a modified version of the [tblite](https://github.com/tblite/tblite) library. Binaries are available for Linux, Windows, and ARM-based macOS.

Extract the tarball or ZIP archive for your operating system from the `binaries` directory into a location on your `PATH` (for example, `$HOME/bin`). No external parameter file is required.

> [!WARNING]
> The macOS and Windows packages include shared libraries (`.dylib` or `.dll`) that must be discoverable at runtime. In most cases, adding the `bin` directory to your `PATH` is sufficient. On macOS, external `.dylib` files may trigger security warnings.

## Usage

The `xtb` binary is largely equivalent to the current [bleeding-edge xtb version](https://github.com/grimme-lab/xtb/releases/tag/bleed). To run a g-xTB calculation, simply add the `--gxtb` flag. For general `xtb` usage, see the [documentation](https://xtb-docs.readthedocs.io/en/latest/) or run `xtb --help`. 

### Singlepoint calculation

```bash
xtb struc.xyz --gxtb
```

You can specify the total charge with `--chrg` or a `.CHRG` file (default: neutral), and the number of unpaired electrons with `--uhf` or a `.UHF` file (default: singlet for systems with an even number of electrons, doublet for systems with an odd number of electrons).

Open-shell g-xTB calculations use an unrestricted wavefunction, which is always triggered if `--uhf` or a `.UHF` file is present, even if the number of unpaired electrons is set to zero. 

Control the verbosity of printed properties with `--verbose` and `--silent`.

### Geometry optimization

```bash
xtb struc.xyz --gxtb --opt
```

Constraints, scans, molecular dynamics, and related features can be controlled through the `xcontrol` input file. To print the gradient, add the `--grad` flag. 

### Numerical Hessian

```bash
xtb struc.xyz --gxtb --hess
```

Computes a numerical Hessian based on the analytic gradient. For the numerical derivative, tight SCF convergence is beneficial. `--acc` defines a multiplicative factor for the convergence criteria (recommended for hessians: 0.1–0.01).

> [!WARNING]
> The macOS binary has a problem with the diagonalization during parallel numerical Hessian computations (serial runs are not affected). 

### Molden file

```bash
xtb struc.xyz --gxtb --molden
```

Writes a Cartesian Molden file including basis-set and orbital information for visualization and post-processing.

### Solvation

g-xTB currently has only limited solvation support through either the generalized Born model with finite epsilon (`--gbe`) or the domain-decomposition conductor-like screening model (`--cosmo`). Both models include only the electrostatic contribution to the solvation free energy. A properly parameterized solvation model for g-xTB is already under development.

```bash
xtb struc.xyz --gxtb --gbe toluene
xtb struc.xyz --gxtb --cosmo water
```

> [!WARNING]
> Preliminary testing indicates that the current form of `gbe` can be unstable during geometry optimization. In addition, the current `cosmo` implementation does not include the solute response to the solvent polarization, which makes the gradient inconsistent. Both models should therefore be used with caution, and results should be evaluated carefully.

## ⚙️ Known limitations

At present, not all `xtb` features are supported for g-xTB. Known limitations include:

- the aISS docking module (`dock`)
- orbital localization (`--lmo`) 
- point charge embedding (`$pcem`)
- cube file generation (`$cube`)

If you encounter another limitation or a bug, please open an issue!
