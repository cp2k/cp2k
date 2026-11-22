# Band Structure Plotting and Export Tool

`plot_band_structure.py` parses CP2K `.bs` files produced by `&DFT%PRINT%BAND_STRUCTURE` (both
standard Kohn-Sham band structures and `&UNFOLDING` output, including spin-polarized runs, multiple
`&KPOINT_SET` blocks, and `ATOMS` groups), renders publication-ready plots, and exports tidy CSV
tables for custom post-processing.

## Requirements

- Python 3.9+
- `numpy`
- `matplotlib` (only required when generating images; `--csv` works with `numpy` alone)

## Common Usage

Plot a standard Kohn-Sham band structure referenced to the Fermi level from `cp2k.out`:

```bash
python3 plot_band_structure.py bands.bs --cp2k-out cp2k.out --emin -3 --emax 3 -o bands.png
```

Plot an unfolded spectral function $A(\mathbf{k}, E)$ (default `--mode spectrum`):

```bash
python3 plot_band_structure.py unfold.bs --cp2k-out cp2k.out --emin -3 --emax 3 --sigma 0.03 -o unfold.png
```

Plot a single `ATOMS` group (for example, dopant contribution in group 2):

```bash
python3 plot_band_structure.py unfold.bs --cp2k-out cp2k.out --spin 1 --group 2 --emin -3 --emax 3 -o dopant.png
```

Plot discrete colored band lines with sparse dots (useful for atom-projected supercell bands with
`TRANSFORMATION_MATRIX 1 0 0 0 1 0 0 0 1`):

```bash
python3 plot_band_structure.py supercell.bs --cp2k-out cp2k.out \
  --mode lines --spin 1 --group 2 --vmax 0.25 --cmap RdYlBu_r \
  --emin -2.5 --emax 2.5 -o supercell_lowdin.png
```

Export parsed band data (including `k_dist`, shifted energies, and optional `weight` / `w_grp*`
columns) to CSV without plotting:

```bash
python3 plot_band_structure.py unfold.bs --cp2k-out cp2k.out --csv unfold.csv
```

## Running Tests

```bash
python3 -m unittest discover -s tools/band_structure -v
```
