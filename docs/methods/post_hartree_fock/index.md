# Post Hartree-Fock

```{toctree}
---
titlesonly:
maxdepth: 2
---
preliminaries
mp2
rpa
low-scaling
```

Post-Hartree-Fock methods in CP2K add wavefunction-based correlation, quasiparticle, or response
corrections on top of a converged reference calculation. The reference is usually Hartree-Fock, a
hybrid functional, or semilocal DFT, depending on the target method and property.

Start with the [preliminaries](preliminaries) page when setting up MP2, RPA, or related methods. It
explains the common choices for primary and RI basis sets, pseudopotentials, reference orbitals,
memory layout, and GPW-based integral grids. The method-specific pages then cover canonical and
RI-MP2, RPA and SOS-MP2, and the low-scaling implementations for larger systems.

These methods are considerably more expensive than semilocal DFT. For production work, converge the
reference calculation first, then check basis-set size, auxiliary basis or auto-generated RI basis
settings, quadrature parameters, and memory distribution.

```{youtube} wyux20qVlck
---
align: center
privacy_mode:
---
```
