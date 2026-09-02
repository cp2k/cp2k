/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*
 * C interface to libsymmetrix (github.com/wcwitt/symmetrix), used by
 * symmetrix_api.F to drive MACE potentials through symmetrix's bespoke kernels.
 *
 * symmetrix has no torch dependency, so this file is compiled in-tree with the
 * regular CP2K C++ compiler and linked against libsymmetrix.a (+ sphericart +
 * BLAS) -- no ABI shim is needed.
 *
 * Contract (see symmetrix `MACE::compute_node_energies_forces`): the engine
 * supplies a full neighbour list in CSR form (num_neigh per central atom, then
 * flattened neigh_indices/neigh_types/displacement-vectors/distances).
 * symmetrix fills per-atom node_energies and PER-EDGE node_forces; the total
 * energy is the sum of node_energies and each edge force is scattered +to the
 * neighbour and -from the central atom. All quantities are in the model's
 * native units (typically Angstrom / eV) -- unit conversion is the Fortran
 * caller's job.
 */

#if defined(__SYMMETRIX)

#include <cstdlib>
#include <cstring>
#include <span>
#include <string>
#include <vector>

// libsymmetrix exposes two interchangeable MACE evaluators with the same load
// /cutoff/types/compute surface: a host-only `MACE` (std::span, std::vector)
// and a Kokkos `MACEKokkos` (device Kokkos::View) selected when libsymmetrix is
// built with SYMMETRIX_KOKKOS=ON. Build CP2K with __SYMMETRIX_KOKKOS to link
// the GPU evaluator; all neighbour-list assembly, force scatter and virial
// below are shared host code -- only the model type and the host<->device
// staging differ.
#if defined(__SYMMETRIX_KOKKOS)
#include "mace_kokkos.hpp"
#include "offload/offload_library.h"
// upstream symmetrix made MACEKokkos a template<Precision>; models are float64.
using SymmetrixModel = MACEKokkos<double>;
#else
#include "mace.hpp"
using SymmetrixModel = MACE;
#endif

static void symmetrix_set_error(int *status, char *err_msg, int err_len,
                                const std::string &msg) {
  *status = 1;
  if (err_msg != nullptr && err_len > 0) {
    std::strncpy(err_msg, msg.c_str(), err_len - 1);
    err_msg[err_len - 1] = '\0';
  }
}

#if defined(__SYMMETRIX_KOKKOS)
// Kokkos must be live before any View is allocated or MACEKokkos is built. Init
// once, lazily, binding this rank to its CP2K offload device assignment (one
// GPU per rank). Kokkos cannot be re-initialized after finalize while CP2K
// force environments come and go, so teardown has to wait for process exit;
// the atexit handler runs after MPI_Finalize, which Kokkos tolerates as long
// as no Kokkos calls are made in between.
static void symmetrix_ensure_kokkos() {
  if (!Kokkos::is_initialized()) {
    Kokkos::InitializationSettings settings;
    const int device_id = offload_get_chosen_device();
    if (device_id >= 0) {
      settings.set_device_id(device_id);
    }
    Kokkos::initialize(settings);
    std::atexit([]() {
      if (Kokkos::is_initialized()) {
        Kokkos::finalize();
      }
    });
  }
}

// Copy a host array into a freshly allocated device View. The incoming pointer
// is wrapped in an unmanaged host view so deep_copy is one contiguous transfer.
template <typename T>
static Kokkos::View<T *> symmetrix_upload(const char *label, const T *host,
                                          int n) {
  Kokkos::View<T *> dev(label, n);
  Kokkos::View<const T *, Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      host_view(host, n);
  Kokkos::deep_copy(dev, host_view);
  return dev;
}

// Copy a device View back into a host vector (one contiguous transfer).
static void symmetrix_download(const Kokkos::View<double *> &dev,
                               std::vector<double> &host) {
  host.resize(dev.extent(0));
  Kokkos::View<double *, Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      host_view(host.data(), host.size());
  Kokkos::deep_copy(host_view, dev);
}
#endif

extern "C" {

/*******************************************************************************
 * \brief Load a symmetrix MACE model from a .json file.
 ******************************************************************************/
void symmetrix_c_model_load(const char *path, void **model_out, int *status,
                            char *err_msg, int err_len) {
  *status = 0;
  *model_out = nullptr;
  try {
#if defined(__SYMMETRIX_KOKKOS)
    symmetrix_ensure_kokkos();
#endif
    *model_out = new SymmetrixModel(std::string(path));
  } catch (const std::exception &e) {
    symmetrix_set_error(status, err_msg, err_len,
                        std::string("symmetrix: failed to load model '") +
                            path + "': " + e.what());
  } catch (...) {
    symmetrix_set_error(status, err_msg, err_len,
                        std::string("symmetrix: failed to load model '") +
                            path + "': unknown C++ exception");
  }
}

/*******************************************************************************
 * \brief Neighbour cutoff (model native length unit, typically Angstrom).
 *
 * Plain field read; cannot throw.
 ******************************************************************************/
double symmetrix_c_model_cutoff(void *model) {
  return static_cast<SymmetrixModel *>(model)->r_cut;
}

/*******************************************************************************
 * \brief Number of atomic types (elements) supported by the model.
 *
 * Plain size query; cannot throw.
 ******************************************************************************/
int symmetrix_c_model_num_types(void *model) {
  return static_cast<int>(
      static_cast<SymmetrixModel *>(model)->atomic_numbers.size());
}

/*******************************************************************************
 * \brief Copy the supported atomic numbers (Z) into out_types.
 *
 * The Kokkos staging copy can throw, so this carries the status protocol.
 ******************************************************************************/
void symmetrix_c_model_types(void *model, int32_t *out_types, int *status,
                             char *err_msg, int err_len) {
  *status = 0;
  try {
    const auto &an = static_cast<SymmetrixModel *>(model)->atomic_numbers;
#if defined(__SYMMETRIX_KOKKOS)
    // atomic_numbers is a device Kokkos::View; stage it through a host mirror.
    auto mirror = Kokkos::create_mirror_view(an);
    Kokkos::deep_copy(mirror, an);
    for (size_t i = 0; i < mirror.extent(0); i++) {
      out_types[i] = static_cast<int32_t>(mirror(i));
    }
#else
    for (size_t i = 0; i < an.size(); i++) {
      out_types[i] = static_cast<int32_t>(an[i]);
    }
#endif
  } catch (const std::exception &e) {
    symmetrix_set_error(status, err_msg, err_len,
                        std::string("symmetrix: atomic-types query failed: ") +
                            e.what());
  } catch (...) {
    symmetrix_set_error(status, err_msg, err_len,
                        "symmetrix: atomic-types query failed: unknown C++ "
                        "exception");
  }
}

/*******************************************************************************
 * \brief Evaluate total energy and per-atom forces, over an atom-decomposed
 *        node set.
 *
 * Nodes are ordered [ owned | ghost | outer ]:
 *   [0, n_local)        atoms this rank owns   -> layer 1 + readouts
 *   [0, n_layer0)       owned + 1-hop ghosts   -> layer 0
 *   [0, n_total)        + 2-hop "outer" atoms  -> array extent only
 * A ghost needs its own complete edge list for a correct H1, and those edges
 * reach two hops out; the outer atoms are therefore edge TARGETS only - they
 * carry positions/types and receive forces, but are never nodes and never have
 * features computed. (LAMMPS pays for these with a 2*r_cut communication
 * cutoff; here they are free, because FIST replicates every coordinate.)
 * Safe because neigh_indices is only dereferenced inside compute_Phi1(n_local),
 * where every target is guaranteed < n_layer0.
 *
 * MACE is a 2-layer message-passing model, so an exact site energy for an owned
 * atom i needs the *layer-1 features* of every atom within r_cut of i, and each
 * of those needs its own complete edge list. Hence the split, exactly as in
 * symmetrix's own LAMMPS pair style (pair_symmetrix_mace.cpp,
 * compute_no_mpi_message_passing):
 *
 *   layer 0 (R0/A0/A0_scaled/M0/H1)      -> n_total  (owned + ghosts)
 *   layer 1 (R1/Phi1/A1/A1_scaled/M1/H2) -> n_local  (owned only)
 *   readouts + ZBL                       -> n_local  (owned only)
 *
 * The reverse pass mirrors it. The `false` accumulate flags on reverse_H2 /
 * reverse_A1_scaled / reverse_Phi1 are load-bearing: they make the two layers'
 * adjoint contributions accumulate rather than clobber. Do not "clean up".
 *
 * ENERGY is masked to owned nodes (compute_readouts(n_local) only ever fills
 * node_energies[0, n_local)), so ghosts cannot double-count.
 * FORCES and the VIRIAL are deliberately NOT masked: every edge of every node,
 * ghost nodes included, contributes. The per-rank result is then exactly
 * d(sum_{i in owned} E_i)/dr_j, a genuine partial derivative; CP2K's existing
 * cross-rank sum in fist_force.F assembles the exact total. This is why the
 * caller must NOT divide by num_pe, unlike the TorchScript MACE/NequIP paths in
 * manybody_e3nn.F which evaluate the whole system on every rank.
 *
 * With n_local == n_total this reduces to exactly the sequence inside
 * libsymmetrix's own compute_node_energies_forces(), i.e. it is bit-identical
 * to the old monolithic call. That equivalence is the P=1 regression gate.
 *
 * Neighbour list (CSR, grouped by central atom, owned nodes' edges first):
 *   node_types[n_total]      symmetrix 0-based type index per atom
 *   num_neigh[n_total]       number of edges whose central atom is this atom
 *   neigh_indices[n_edges]   neighbour (j) packed atom index per edge
 *   neigh_types[n_edges]     symmetrix type index of the neighbour
 *   edge_centre[n_edges]     central (i) packed atom index per edge
 *   xyz[3*n_edges]           displacement vector r_j - r_i per edge
 *   r[n_edges]               |displacement| per edge
 * out_forces[3*n_total]; out_virial filled when do_virial is nonzero.
 * Native model units throughout.
 ******************************************************************************/
void symmetrix_c_eval_mp(void *model, int n_local, int n_layer0, int n_total,
                         const int32_t *node_types, int n_edges,
                         const int32_t *num_neigh, const int32_t *neigh_indices,
                         const int32_t *neigh_types, const int32_t *edge_centre,
                         const double *xyz, const double *r, int do_virial,
                         double *out_energy, double *out_forces,
                         double *out_virial, int *status, char *err_msg,
                         int err_len) {
  *status = 0;
  try {
    auto *mace = static_cast<SymmetrixModel *>(model);

    // symmetrix takes int (32-bit on all CP2K targets); copy into int vectors
    // to be ABI-safe and explicit, and to stage the Kokkos device upload below.
    std::vector<int> nt(node_types, node_types + n_total);
    std::vector<int> nn(num_neigh, num_neigh + n_total);
    std::vector<int> ni(neigh_indices, neigh_indices + n_edges);
    std::vector<int> ntyp(neigh_types, neigh_types + n_edges);

    // node_energies[n_local] and per-edge node_forces[3*n_edges], read back on
    // the host so the scatter/virial below is identical for both evaluators.
    const double *node_energies = nullptr;
    const double *node_forces = nullptr;
#if defined(__SYMMETRIX_KOKKOS)
    std::vector<double> ne_host, nf_host;
    {
      // One upload / one download around the WHOLE layer sequence - never
      // round-trip per layer, that would dwarf the compute.
      auto d_nt = symmetrix_upload<int>("node_types", nt.data(), n_total);
      auto d_nn = symmetrix_upload<int>("num_neigh", nn.data(), n_total);
      auto d_ni = symmetrix_upload<int>("neigh_indices", ni.data(), n_edges);
      auto d_ntyp = symmetrix_upload<int>("neigh_types", ntyp.data(), n_edges);
      auto d_xyz = symmetrix_upload<double>("xyz", xyz, 3 * n_edges);
      auto d_r = symmetrix_upload<double>("r", r, n_edges);

      if (mace->node_energies.size() < static_cast<size_t>(n_local))
        Kokkos::realloc(mace->node_energies, n_local);
      Kokkos::deep_copy(mace->node_energies, 0.0);
      if (mace->node_forces.size() < static_cast<size_t>(3 * n_edges))
        Kokkos::realloc(mace->node_forces, 3 * n_edges);
      Kokkos::deep_copy(mace->node_forces, 0.0);

      if (mace->has_zbl)
        mace->zbl.compute_ZBL(n_local, d_nt, d_nn, d_ntyp, mace->atomic_numbers,
                              d_r, d_xyz, mace->node_energies,
                              mace->node_forces);

      mace->compute_Y(d_xyz);

      mace->compute_R0(n_layer0, d_nt, d_nn, d_ntyp, d_r);
      mace->compute_A0(n_layer0, d_nt, d_nn, d_ntyp);
      mace->compute_A0_scaled(n_layer0, d_nt, d_nn, d_ntyp, d_r);
      mace->compute_M0(n_layer0, d_nt);
      mace->compute_H1(n_layer0);

      mace->compute_R1(n_local, d_nt, d_nn, d_ntyp, d_r);
      mace->compute_Phi1(n_local, d_nn, d_ni);
      mace->compute_A1(n_local);
      mace->compute_A1_scaled(n_local, d_nt, d_nn, d_ntyp, d_r);
      mace->compute_M1(n_local, d_nt);
      mace->compute_H2(n_local, d_nt);

      mace->compute_readouts(n_local, d_nt);

      mace->reverse_H2(n_local, d_nt, false);
      mace->reverse_M1(n_local, d_nt);
      mace->reverse_A1_scaled(n_local, d_nt, d_nn, d_ntyp, d_xyz, d_r);
      mace->reverse_A1(n_local);
      mace->reverse_Phi1(n_local, d_nn, d_ni, d_xyz, d_r, false, false);

      mace->reverse_H1(n_layer0);
      mace->reverse_M0(n_layer0, d_nt);
      mace->reverse_A0_scaled(n_layer0, d_nt, d_nn, d_ntyp, d_xyz, d_r);
      mace->reverse_A0(n_layer0, d_nt, d_nn, d_ntyp, d_xyz, d_r);

      symmetrix_download(mace->node_energies, ne_host);
      symmetrix_download(mace->node_forces, nf_host);
    }
    node_energies = ne_host.data();
    node_forces = nf_host.data();
#else
    const auto s_nt = std::span<const int>(nt);
    const auto s_nn = std::span<const int>(nn);
    const auto s_ni = std::span<const int>(ni);
    const auto s_ntyp = std::span<const int>(ntyp);
    const auto s_xyz = std::span<const double>(xyz, 3 * n_edges);
    const auto s_r = std::span<const double>(r, n_edges);

    mace->node_energies.assign(n_local, 0.0);
    mace->node_forces.assign(3 * n_edges, 0.0);

    if (mace->has_zbl)
      mace->zbl.compute_ZBL(n_local, s_nt, s_nn, s_ntyp, mace->atomic_numbers,
                            s_r, s_xyz, mace->node_energies, mace->node_forces);

    mace->compute_Y(s_xyz);

    mace->compute_R0(n_layer0, s_nt, s_nn, s_ntyp, s_r);
    mace->compute_A0(n_layer0, s_nt, s_nn, s_ntyp);
    mace->compute_A0_scaled(n_layer0, s_nt, s_nn, s_ntyp, s_r);
    mace->compute_M0(n_layer0, s_nt);
    mace->compute_H1(n_layer0);

    mace->compute_R1(n_local, s_nt, s_nn, s_ntyp, s_r);
    mace->compute_Phi1(n_local, s_nn, s_ni);
    mace->compute_A1(n_local);
    mace->compute_A1_scaled(n_local, s_nt, s_nn, s_ntyp, s_r);
    mace->compute_M1(n_local, s_nt);
    mace->compute_H2(n_local, s_nt);

    mace->compute_readouts(n_local, s_nt);

    mace->reverse_H2(n_local, s_nt, false);
    mace->reverse_M1(n_local, s_nt);
    mace->reverse_A1_scaled(n_local, s_nt, s_nn, s_ntyp, s_xyz, s_r, false);
    mace->reverse_A1(n_local);
    mace->reverse_Phi1(n_local, s_nn, s_ni, s_xyz, s_r, false, false);

    mace->reverse_H1(n_layer0);
    mace->reverse_M0(n_layer0, s_nt);
    mace->reverse_A0_scaled(n_layer0, s_nt, s_nn, s_ntyp, s_xyz, s_r);
    mace->reverse_A0(n_layer0, s_nt, s_nn, s_ntyp, s_xyz, s_r);

    node_energies = mace->node_energies.data();
    node_forces = mace->node_forces.data();
#endif

    // energy: owned nodes only (readouts never filled the ghost entries)
    double energy = 0.0;
    for (int a = 0; a < n_local; a++) {
      energy += node_energies[a];
    }
    *out_energy = energy;

    // per-edge node_forces scattered: +to neighbour j, -from central atom i.
    // Every edge contributes, ghost-node edges included - see the note above.
    std::memset(out_forces, 0, sizeof(double) * 3 * n_total);
    for (int e = 0; e < n_edges; e++) {
      const int i = edge_centre[e];
      const int j = neigh_indices[e];
      for (int d = 0; d < 3; d++) {
        const double f = node_forces[3 * e + d];
        out_forces[3 * j + d] += f;
        out_forces[3 * i + d] -= f;
      }
    }

    // virial = + sum_edges node_force (outer) displacement, matching CP2K's
    // pv_nonbond convention (energy units). Summed over ALL edges, like the
    // forces. NOTE: with an atom-decomposed node set the per-rank virial is
    // NOT symmetric - each rank holds a lopsided subset of the edges. Only the
    // cross-rank sum (fist_force.F) is. Do not "fix" the asymmetry here.
    std::memset(out_virial, 0, sizeof(double) * 9);
    if (do_virial != 0) {
      double vir[9] = {0.0};
      for (int e = 0; e < n_edges; e++) {
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            vir[3 * a + b] += node_forces[3 * e + a] * xyz[3 * e + b];
          }
        }
      }
      std::memcpy(out_virial, vir, sizeof(double) * 9);
    }
  } catch (const std::exception &e) {
    symmetrix_set_error(status, err_msg, err_len,
                        std::string("symmetrix: evaluation failed: ") +
                            e.what());
  } catch (...) {
    symmetrix_set_error(status, err_msg, err_len,
                        "symmetrix: evaluation failed: unknown C++ exception");
  }
}

/*******************************************************************************
 * \brief Free a model handle.
 ******************************************************************************/
void symmetrix_c_model_release(void *model) {
  delete static_cast<SymmetrixModel *>(model);
}

/*******************************************************************************
 * \brief Whether this binary evaluates symmetrix on a Kokkos device backend
 *        (build-time choice, CP2K_SYMMETRIX_KOKKOS) rather than on the host.
 ******************************************************************************/
int symmetrix_c_is_kokkos(void) {
#if defined(__SYMMETRIX_KOKKOS)
  return 1;
#else
  return 0;
#endif
}

} // extern "C"

#endif // __SYMMETRIX
