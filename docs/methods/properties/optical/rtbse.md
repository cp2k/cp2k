# Real Time Bethe-Salpeter Propagation

Instead of using TDDFT functionals, an approximation to the self-energy is employed to calculate the
time dependent behaviour of the density matrix \[[Attaccalite2011](http://dx.doi.org/10.1103/PhysRevB.84.245110)\].
This requires a previous determination of the screened Coulomb potential, done via the bandstructure 
[GW](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW) calculation.

The single particle density matrix $\hat{\rho}$ is propagated following the equation of motion

$$ \frac{\mathrm{d} \hat{\rho}}{\mathrm{d} t} = -i [\hat{H}(t), \hat{\rho}(t)] $$

This equation is solved by steps

$$ \hat{\rho} (t + \Delta t) = \mathrm{e} ^ {- i \hat{H} (t+\Delta t) \Delta t/2} \mathrm{e} ^ {-i \hat{H}(t) \Delta t/2}
\hat{\rho (t)} \mathrm{e} ^ {i \hat{H}(t) \Delta t/2} \mathrm{e} ^ {i \hat{H} (t + \Delta t) \Delta t/2}$$

which is called the _enforced time reversal scheme_\[[Castro2004](https://doi.org/10.1063/1.1774980)\].
The effective Hamiltonian is given as

$$ \hat{H}(t) = \hat{h}^{G0W0} + \hat{U} (t) +
\hat{V}^{\mathrm{Hartree}} [\hat{\rho}(t)] - \hat{V}^{\mathrm{Hartree}} [\hat{\rho}_0] +
\hat{\Sigma}^{\mathrm{COHSEX}}[\hat{\rho}(t)] - \hat{\Sigma}^{\mathrm{COHSEX}}[\hat{\rho}_0]
$$

where $\hat{\rho}_0$ is the density matrix determined from the molecular orbitals
used in [GW](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW) and $\hat{U}(t)$ is the
external applied field.

## Excitation scheme

Without the external field $\hat{U}(t)$, the density matrix only rotates in phase but does not
produce any measurable dynamics. The excitation of the dynamics can be done either by a
real time pulse (i.e. at each point, $\hat{U}(t)$ follows form due to some finite time dependent
field $\vec{E}(t)$) or by an infinitely sharp delta pulse, which we can understand as the limit of
$\vec{E}(t) \to I \vec{e} \delta(t)$, where $I$ is the delta pulse intensity and $\vec{e}$ its direction.

## Observables

The dynamics can be traced through time with electric dipole moment associated with the density matrix

$$ \mu_i(t) = \mathrm{Tr} (\hat{\rho}(t) \hat{x}_i)
$$

The electric polarizability (which is related to the photon absorption spectrum) is then
determined as

$$ \alpha_{ij} (\omega) = \frac{\mu_i(\omega)}{E_j(\omega)}
$$

where we Fourier transformed to the frequency domain. In order to model the Fourier transform
of infinitely oscillating dipole moments, we introduce a damping factor $\gamma$\[[Müller2020](https://doi.org/10.1002/jcc.26412)\]

$$ \mu_i(\omega) = \int _ 0 ^ T \mathrm{d}t \mathrm{e}^{-\gamma t} \mathrm{e} ^ {i \omega t} \mu_i(t)
$$

One can easily verify that for real FT of the applied field, this leads Lorentzian peaks
at the frequencies of the oscillations of the moments present in the imaginary part of the
corresponding polarizability element.

Use [RTP_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.RTP_METHOD) to start the
calculation by setting it to TDAGW.
