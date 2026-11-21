#!/usr/bin/env python3

# author: BSG

import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Run RIXS convolution script.")
    parser.add_argument(
        "--filename", type=str, help="Input data file from cp2k (*.rixs)"
    )
    parser.add_argument("--v_states", type=int, help="Number of valence excited states")
    parser.add_argument(
        "--gamma", type=float, default=0.075, help="Broadening width in eV"
    )
    parser.add_argument(
        "--kappa", type=float, default=60.0, help="Scattering angle in degrees"
    )
    return parser.parse_args()


# hbar = 6.6260689 * 10**-34 / (2 * np.pi)
hbar = 1  # a.u.


def get_sigma_tensor(f_state, w_in, w_0i, w_if, mu_0i, mu_if, Gamma):
    """
    Calculate absorption cross-section tensor of final valence state f:
    sigma_{alpha beta gamma delta}^f(w)

    Parameters:
    - f_state: int, index of the final valence state
    - w_in: float or array, incident photon frequencies
    - w_0i: array (c_states,)
    - w_if: array (c_states, v_states)
    - mu_0i: array (c_states, 3) dipole moments (i->0)
    - mu_if: array (c_states, v_states, 3) dipole moments (i->f)
    - Gamma: float, broadening parameter

    Returns:
    - sigma: array (3,3,3,3), or (n_energies,3,3,3,3) for array input
    """
    energies = np.atleast_1d(w_in)
    prefac = (w_if[:, f_state] * w_0i) ** 2 / (
        (energies[:, None] - w_0i[None, :]) ** 2 + Gamma**2
    )
    transition = mu_if[:, f_state, :]
    moments = np.einsum("ia,ib,ic,id->iabcd", transition, transition, mu_0i, mu_0i)
    sigma = np.einsum("ri,iabcd->rabcd", prefac, moments)
    return sigma[0] if np.ndim(w_in) == 0 else sigma


def parse_spectrum_file(filename, v_states):
    """
    Parameters:
    - filename: str, path to the spectrum file
    - v_states: int, number of valence excited states
    """

    w_0i = []
    mu_0i = []
    w_if = []
    mu_if = []
    w_f0 = []

    with open(filename, "r") as f:
        lines = f.readlines()

    in_core_block = False
    in_valence_block = False

    count = 0

    for line in lines:
        line = line.strip()

        if "Excitation from ground-state" in line:
            in_core_block = True
            in_valence_block = False
            continue
        elif "Emission from core-excited state" in line:
            in_core_block = False
            in_valence_block = True
            continue
        elif line.startswith("=") or line.startswith("w_0i") or not line:
            continue

        tokens = line.split()

        if in_core_block and len(tokens) >= 5:
            e0i, mu_x, mu_y, mu_z = map(float, tokens[0:4])
            w_0i.append(e0i / hbar)
            mu_0i.append([mu_x, mu_y, mu_z])

        elif in_valence_block and len(tokens) >= 6:
            count += 1
            ei, ef = map(float, tokens[0:2])
            mu_x, mu_y, mu_z = map(float, tokens[2:5])
            mu_if.append((mu_x, mu_y, mu_z))
            w_if.append((ef - ei) / hbar)
            if count <= v_states:
                w_f0.append(ef / hbar)

    xas_w_0i = np.array(w_0i)
    xas_mu_0i = np.array(mu_0i)
    rixs_w_if = np.array(w_if)
    tddft_w_f0 = np.array(w_f0)
    c_states = len(mu_0i)

    # sanity check
    expected_total = c_states * v_states
    if len(mu_if) != expected_total:
        raise ValueError(
            f"Expected {expected_total} spectrum lines, but found {len(mu_if)}"
        )
    rixs_w_if = rixs_w_if.reshape((c_states, v_states))

    rixs_mu_if = np.zeros((c_states, v_states, 3))

    for idx, (mx, my, mz) in enumerate(mu_if):
        i = idx // v_states
        f = idx % v_states
        if i < c_states and f < v_states:
            rixs_mu_if[i, f] = [mx, my, mz]

    return xas_w_0i, xas_mu_0i, rixs_w_if, rixs_mu_if, tddft_w_f0


def convolve_sigma_averaged_tensor(
    v_states, w_f0, w_0i, w_if, mu_0i, mu_if, w_in_window, kappa_rad, Gamma
):
    """
    Compute the RIXS spectrum with custom broadening.

    Parameters:
    - v_states: int, number of valence excited states
    - w_f0: (v_states,) energy loss values (valence excitation energies)
    - w_0i: (c_states,) core excitation energies
    - w_if: (c_states, v_states) emission energies
    - mu_0i: (c_states, 3) ground-to-core dipole moments
    - mu_if: (c_states, v_states, 3) core-to-valence dipole moments
    - kappa_rad: float, scattering angle in radians
    - Gamma: float, broadening width

    Returns:
    - 2D array representing the RIXS spectrum
    """
    A = 3 + np.cos(kappa_rad) ** 2
    B = 0.5 * (1 - 3 * np.cos(kappa_rad) ** 2)

    w_out_window = np.linspace(0, np.max(w_f0), len(w_in_window))
    N_incident = len(w_in_window)
    N_loss = len(w_out_window)
    spectrum = np.zeros((N_incident, N_loss))

    for f in range(v_states):
        sigma = get_sigma_tensor(f, w_in_window, w_0i, w_if, mu_0i, mu_if, Gamma)
        delta = Gamma / np.pi / ((w_f0[f] - w_out_window) ** 2 + Gamma**2)
        val = A * np.einsum("raabb->r", sigma)
        val += B * (np.einsum("rabab->r", sigma) + np.einsum("rabba->r", sigma))
        spectrum += val[:, None] * delta[None, :] / 30.0

    return spectrum / np.max(spectrum)


if __name__ == "__main__":
    args = parse_args()
    filename = args.filename
    v_states = args.v_states
    kappa_deg = args.kappa
    kappa_rad = np.radians(kappa_deg)
    Gamma = args.gamma

    w_0i, mu_0i, w_if, mu_if, w_f0 = parse_spectrum_file(filename, v_states)

    w_in = w_0i[0]  # incident photon energy = first xas energy
    e_window = 0.2
    n_points = 50
    w_in_window = np.linspace(w_in - e_window, w_in + e_window, n_points)

    spectrum_2d = convolve_sigma_averaged_tensor(
        v_states, w_f0, w_0i, w_if, mu_0i, mu_if, w_in_window, kappa_rad, Gamma
    )

    project_name = filename.split(".")[0]
    np.savetxt(f"{project_name}_rixs.dat", spectrum_2d, fmt="%.6f")
