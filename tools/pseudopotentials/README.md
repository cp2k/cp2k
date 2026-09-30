# SSSP v2.0 import audit

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

The `tests/QS/regtest-upf` inputs combine two different NC potentials with nonlinear core
corrections and the existing UZH-MOLOPT basis sets. They exercise complete energy/force
calculations. Additional checks must cover active GAPW one-center expansions, finite-difference
forces and stress, spin polarization, and basis/grid convergence against an independent UPF
implementation before quantitative SSSP compatibility can be claimed.

Keep three convergence parameters separate: the Gaussian orbital basis, the representation of
tabulated potential/projector data, and the density quadrature. A converged density grid or
agreement of analytic forces with finite differences does not establish the accuracy of the
transformed Hamiltonian. Existing UZH-MOLOPT bases can be used as a starting point, with their
valence partition and convergence checked for the selected UPF dataset.

USPP and PAW are not supported by this NC transformation. Their augmentation functions, overlap
operators, and potential-specific one-center data must be retained and implemented explicitly.
Relabeling either family as NC is invalid.
