# Optical Spectroscopy

```{toctree}
---
titlesonly:
maxdepth: 1
---
tddft
bethe-salpeter
rtbse
```

Optical spectroscopy is a technique used to study the interaction between light and matter. It
involves measuring the absorption, emission, or scattering of light by molecules, atoms, or
materials. The resulting spectra provide valuable information about the electronic structure, energy
levels, and dynamics of the system under investigation.

In optical spectroscopy, one of the key quantities of interest is the excitation energy
($\Omega_n$). This is the energy required to excite a molecule from its ground state to an excited
state. These excitation energies are directly related to the positions and intensities of spectral
lines observed in absorption and emission spectra.

The excitation energies $\Omega_n$ can be computed using various theoretical approaches. Two
commonly used methods are linear-response Time-Dependent Density Functional Theory (LR-TDDFT) and
the linear-response *GW*/Bethe-Salpeter Equation (*GW*/BSE) approach. Both methods can be formulated
using Casida's equation

$$\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad .$$

We abbreviate $A$ and $B$ as matrices with index $A_{ia,jb}$, i.e. they have
$N_\mathrm{occ}N_\mathrm{empty}$ rows and $N_\mathrm{occ}N_\mathrm{empty}$ columns. The matrices $A$
and $B$ are different in TDDFT and *GW*/BSE; for TDDFT they read (for singlet excitations, details
on the TDDFT page)

$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a^\text{DFT}-\varepsilon_i^\text{DFT})\delta_{ij}\delta_{ab} + 
    2v_{ia,jb} + \langle ia|f_\text{xc}(\Omega^{(n)})|jb\rangle \quad ,\\[0.5em]
    B_{ia,jb} &= 2v_{ia,bj} +  \langle ia|f_\text{xc}(\Omega^{(n)})|jb\rangle \quad ,
\end{align}$$

and for *GW*/BSE (details on the *GW*/BSE page):

$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a^{GW}-\varepsilon_i^{GW})\delta_{ij}\delta_{ab} + 
    2v_{ia,jb} - W_{ij,ab} \quad ,\\[0.5em]
    B_{ia,jb} &= 2 v_{ia,bj} - W_{ib,aj} \quad .
\end{align}$$

TDDFT with the common Adiabatic Local Density Approximation (ALDA) or with a hybrid functional (i.e.
PBE0) can be a good choice for calculating excitation energies of molecules. Exceptions include
charge-transfer excitations where the excited electron is transferred over a significant distance
within the molecule. In such cases, range-separated hybrid functionals might be needed.

For solids, the applicability of TDDFT can depend on whether the solid is metallic or has a finite
bandgap. For metals, ALDA often yields good excitation energies. However, for semiconductors and
insulators, ALDA fails because the ALDA xc kernel does not adequately include the Coulomb
interaction between the electron and the hole of the electron-hole pair (exciton) that forms upon
excitation. In contrast, the *GW*/BSE approach is well-suited for computing the excitation energies
of excitons in semiconductors and insulators. *GW*/BSE accounts for the attractive interaction
between the electron and hole in the A-matrix via the screened Coulomb interaction $W_{ij,ab}$. This
inclusion is crucial for accurately describing excitonic effects, which are significant in materials
with a finite bandgap.

Thus, TDDFT with ALDA/hybrid functionals is convenient and computationally less demanding than
*GW*/BSE for molecular systems and metals, *GW*/BSE can describe excitonic effects in semiconductors
and insulators. For a more detailed discussion on TDDFT and *GW*/BSE, we recommend for example C. A.
Ullrich, *Time-Dependent Density-Functional Theory - Concepts and Applications*.
