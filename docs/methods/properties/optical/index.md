# Optical Spectroscopy

```{toctree}
---
titlesonly:
maxdepth: 1
---
tddft
bethe-salpeter
```
Optical spectroscopy is a technique used to study the interaction between light and matter.
It involves measuring the absorption, emission, or scattering of light by molecules, atoms, 
or materials. The resulting spectra provide valuable information about the electronic structure, 
energy levels, and dynamics of the system under investigation.

In optical spectroscopy, one of the key quantities of interest is the excitation energy ($\Omega_n$). 
This is the energy required to excite a molecule from its ground state to an excited state. 
These excitation energies are directly related to the positions and 
intensities of spectral lines observed in absorption and emission spectra.

The excitation energies $\Omega_n$ can be computed using various theoretical approaches. 
Two commonly used methods are Time-Dependent Density Functional Theory (TDDFT) and 
the *GW*/Bethe-Salpeter Equation (*GW*/BSE) approach. Both methods can be formulated using
Casida's equation

$$\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad .$$

We abbreviate $A$ and $B$ as matrices with index $A_{ia,jb}$, i.e. they have
$N_\mathrm{occ}N_\mathrm{empty}$ rows and $N_\mathrm{occ}N_\mathrm{empty}$ columns.
The matrices $A$ and $B$ are different in TDDFT and *GW*/BSE
The entries of
$A$ and

$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a^{GW}-\varepsilon_i^{GW})\delta_{ij}\delta_{ab} + \alpha^\mathrm{S/T}
    v_{ia,jb} - W_{ij,ab}(\omega=0) \quad ,\\
    B_{ia,jb} &= \alpha^\mathrm{(S/T)} v_{ia,bj} - W_{ib,aj}(\omega=0) \quad .
\end{align}$$
