# SSSP v2.0 import audit

For the tested scope, reproducibility instructions, independent references, and remaining
limitations, see [NC_VALIDATION.md](NC_VALIDATION.md).

The pinned archives in `sssp_v2_archives.json` contain the SSSP v2.0 PBE and PBEsol Efficiency and
Precision libraries. The release version is distinct from the UPF file format version (2.0.1).

Fetch, verify, and inventory the four libraries:

```sh
python3 tools/pseudopotentials/sssp_v2_audit.py /path/to/sssp-v2-data
```

Check every NC entry with a built CP2K executable:

```sh
cmake --build /path/to/build --target atom_upf_unittest
python3 tools/pseudopotentials/sssp_v2_audit.py /path/to/sssp-v2-data \
  --check-executable /path/to/build/bin/atom_upf_unittest.psmp
```

The audit checks file hashes, inventory, import, potential lifecycle, and finite Gaussian
transformation data. It does not establish physical accuracy of the Gaussian fits or validate the
SSSP transferability benchmarks. There are 82 NC entries across the four libraries; several occur in
more than one library. The executable also checks analytic radial core-density derivatives against
finite differences, the spherical-harmonic addition theorem through angular momentum five, and
atomic local/nonlocal operators against analytic Gaussian integrals.

## Direct radial projector integrals

Quickstep retains the scalar NC radial projector identities and the full UPF coupling matrix,
preserving the numerical amplitude of the tabulated projectors. The reader removes the radial
prefactor from `PP_BETA` and multiplies `PP_DIJ` by 0.5 to express the complete operator in Hartree.
The unit test checks projector amplitudes separately from the operator: compensating changes to beta
and D alone would conceal a mismatch with future augmentation data. See the
[UPF specification](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
for the radial prefactors and dataset unit conventions. `upf_projector_integrals` integrates
Cartesian polynomials and real spherical harmonics analytically, with a one-dimensional Simpson
quadrature on the supplied radial mesh. Gaussian basis contractions, projector identities, and
off-diagonal radial couplings are preserved. The local potential is also integrated from its
original radial table. Only the NLCC core density retains a Gaussian expansion, fitted with a
column-normalized direct SVD and constrained by zero-density tail samples.

Run the independent analytic integral checks with:

```sh
cmake --build /path/to/build --target upf_projector_integrals_unittest
/path/to/build/bin/upf_projector_integrals_unittest.psmp
```

Synthetic Gaussian radial projectors are compared with the existing analytic Cartesian Gaussian
integral routines through g orbitals, including diffuse primitives, four center separations, linear
and logarithmic meshes, derivatives through second order, and moments through quadrupole order.
These tests verify the integral algebra. They do not measure quadrature convergence for every UPF
dataset or convergence of the molecular orbital basis. Optional UPF filenames on the command line
additionally check native integral finiteness and successive replacement of the radial/angular cache
for each dataset.

## Direct local-potential integrals

`upf_local_integrals` combines the Gaussian product theorem with analytic angular integration and
radial quadrature of the original local potential. It subtracts the Gaussian-core electrostatic
contribution consistently with the core exponent. Every supplied radial sample is retained:
`EPS_PPL` only controls neighbor-list radii, not a truncation of the quadrature mesh. The native
operator is used for Hamiltonian elements, forces, stress, RI integrals, and position moments.

```sh
cmake --build /path/to/build --target upf_local_integrals_unittest
/path/to/build/bin/upf_local_integrals_unittest.psmp
```

The unit test uses independent analytic Gaussian matrix elements, derivatives, and RI references,
with even and odd linear meshes and a logarithmic mesh. A coarse `EPS_PPL` case detects accidental
truncation of the radial integration.

The `tests/QS/regtest-upf` inputs combine two different NC potentials with nonlinear core
corrections and the existing UZH-MOLOPT basis sets. They exercise complete energy/force calculations
in GPW, GAPW_XC, and GAPW, with active one-center corrections in the latter two modes. These compact
regression tests do not establish basis/grid convergence or reproduce the SSSP transferability
benchmarks. Such claims require convergence studies and comparison with an independent UPF
implementation.

Keep three convergence parameters separate: the Gaussian orbital basis, the representation of
tabulated potential/projector data, and the density quadrature. A converged density grid or
agreement of analytic forces with finite differences does not establish the accuracy of the
transformed Hamiltonian. Existing UZH-MOLOPT bases can be used as a starting point, with their
valence partition and convergence checked for the selected UPF dataset.

USPP and PAW are not supported by this NC transformation. Their augmentation functions, overlap
operators, and potential-specific one-center data must be retained and implemented explicitly.
Relabeling either family as NC is invalid.
