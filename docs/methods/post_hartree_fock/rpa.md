# Random-Phase Approximation

```{toctree}
---
titlesonly:
maxdepth: 1
---
rpa-quartic
rpa-kpoints
```

The random-phase approximation (RPA) is a post-Hartree-Fock method for the total energy of
molecules, liquids, surfaces and solids. It starts from a converged DFT calculation and replaces the
exchange-correlation energy of DFT by the exact exchange energy and the RPA correlation energy,
which contains the van der Waals interaction.

## History

RPA was introduced in the early 1950s for the electron gas
([Bohm and Pines 1953](https://doi.org/10.1103/PhysRev.92.609),
[Gell-Mann and Brueckner 1957](https://doi.org/10.1103/PhysRev.106.364)). In solid-state physics, it
has since been the standard description of screening, for example in the dielectric function and in
the screened Coulomb interaction $W$ of the *GW* method. As a total-energy method for molecules and
solids, RPA became feasible around 2000 to 2010
([Furche 2001](https://doi.org/10.1103/PhysRevB.64.195120),
[Harl *et al.* 2010](https://doi.org/10.1103/PhysRevB.81.115126)). Reviews are
[Eshuis *et al.* 2012](https://doi.org/10.1007/s00214-011-1084-8) and
[Ren *et al.* 2012](https://doi.org/10.1007/s10853-012-6570-4).

## RPA total energy

An RPA calculation starts with a converged DFT calculation, commonly with the PBE functional. From
the DFT orbitals and eigenvalues, the RPA total energy is computed once (without self-consistency)
as

$$
E^\text{RPA} = E^\text{DFT} - E^\text{xc} + E^\text{EXX} + E_\text{c}^\text{RPA}\,,
$$

where $E^\text{DFT}$ is the DFT total energy and $E^\text{xc}$ the exchange-correlation energy of
the DFT functional. $E^\text{EXX}$ is the Hartree-Fock exchange energy evaluated with the DFT
orbitals. $E_\text{c}^\text{RPA}$ is the correlation energy from the density response of the DFT
orbitals; it is non-local and contains the van der Waals interaction.

## Why RPA?

Semilocal and hybrid DFT functionals do not contain the van der Waals interaction; it is usually
added as a correction with empirical parameters or as a non-local functional. RPA contains the van
der Waals interaction without any parameter and on the same footing as covalent, ionic and metallic
bonding. RPA is therefore well suited for systems with several types of bonding:

- **Adsorption on surfaces:** For CO on transition-metal surfaces, RPA gives accurate surface
  energies and adsorption energies at the same time, which semilocal functionals do not
  ([Schimka *et al.* 2010](https://doi.org/10.1038/nmat2806)). For graphene on Ni(111), RPA
  describes the competition between chemical and van der Waals bonding
  ([Mittendorfer *et al.* 2011](https://doi.org/10.1103/PhysRevB.84.201401)).
- **Layered materials:** RPA reproduces the interlayer binding of graphite
  ([Lebègue *et al.* 2010](https://doi.org/10.1103/PhysRevLett.105.196401)).
- **Liquids:** Monte Carlo simulations with RPA total energies in CP2K give a density of liquid
  water close to experiment ([Del Ben *et al.* 2015](https://doi.org/10.1063/1.4927325)).

For lattice constants and bulk moduli of bulk solids, RPA is more accurate than PBE
([Harl *et al.* 2010](https://doi.org/10.1103/PhysRevB.81.115126)), but meta-GGA functionals like
SCAN reach a similar accuracy at much lower cost
([Sun *et al.* 2016](https://doi.org/10.1038/nchem.2535)). For bulk solids alone, RPA is thus not
necessarily needed.

**Comparison to MP2 and coupled cluster.** Second-order Møller-Plesset perturbation theory (MP2)
contains no screening of the Coulomb interaction. The MP2 correlation energy diverges for metals,
and MP2 overestimates the binding of solids with a small band gap
([Grüneis *et al.* 2010](https://doi.org/10.1063/1.3466765)). RPA includes screening and can be used
for insulators, semiconductors and metals. Coupled-cluster theory, CCSD(T), is more accurate than
RPA for insulators ([Booth *et al.* 2013](https://doi.org/10.1038/nature11770)), but much more
expensive.

**Computation time.** An RPA calculation is more expensive than a DFT calculation with a semilocal
functional. As a rough guide, the computation time is similar to or larger than that of a
hybrid-functional calculation.

**Limitations.** Binding energies and atomization energies from RPA are typically too small
([Harl *et al.* 2010](https://doi.org/10.1103/PhysRevB.81.115126)). Corrections like renormalized
single excitations ([Ren *et al.* 2011](https://doi.org/10.1103/PhysRevLett.106.153003)) and
second-order screened exchange \[[](#Grueneis2009)\] reduce this error. The RPA total energy depends
on the DFT functional of the preceding SCF.

## RPA implementations in CP2K

CP2K contains three RPA implementations:

1. **[RPA for molecules and large supercells](rpa-quartic.md).** $\Gamma$-point only, for example
   for liquid water \[[](#DelBen2013), [](#DelBen2015)\]. Analytical gradients are available
   \[[](#Stein2024)\], as well as corrections beyond RPA (RSE, AXK, SOSEX).
1. **[Low-scaling RPA](low-scaling.md).** $\Gamma$-point only, cubic scaling, for large systems
   \[[](#Wilhelm2016b), [](#Bussy2023)\].
1. **[RPA with *k*-points (RI-RS)](rpa-kpoints.md).** For crystals with a small unit cell and
   *k*-point sampling. It is based on the real-space resolution of the identity (RI-RS) of Duchemin
   and Blase \[[](#Duchemin2019)\], see also
   [Delesma *et al.* 2024](https://doi.org/10.1063/5.0184406). Nuclear gradients and the stress
   tensor are in preparation.

Double-hybrid functionals contain a fraction of the RPA correlation energy. This rescaling is set
with the keyword [SCALE_RPA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.SCALE_RPA) of the
implementations 1 and 2. The closely related Laplace-transformed scaled-opposite-spin MP2 method is
described in [](sos-mp2.md).
