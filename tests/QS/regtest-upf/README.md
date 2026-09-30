# SSSP v2.0 NC-UPF regression inputs

These inputs use the unmodified PBE Efficiency v2.0 S and Br NC potentials
and the existing UZH-MOLOPT Gaussian basis sets. Both datasets contain a
nonlinear core correction. The SBr2 inputs exercise GPW, GAPW_XC, and GAPW,
including successive UPF kinds, atomic initial densities, and forces.

The source archives and SHA-256 hashes are pinned in
`tools/pseudopotentials/sssp_v2_archives.json`. The original generator and
license notices are retained inside each UPF file.

These regression cases detect implementation failures. Their energies do
not establish convergence to a plane-wave reference or SSSP transferability.
