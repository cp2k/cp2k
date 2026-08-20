/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include <stdio.h>

#include "nnp_gpu_internal.h"
#include "nnp_gpu_state.h"

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

/*
 * Committee neural-network evaluation on the device. The topology, activation
 * functions, atomic energies and packed weights and biases are uploaded once
 * per run (nnp_gpu_upload_nn_packed_c); the per-step launcher
 * (nnp_gpu_nn_atomwise_with_G_dev_c) evaluates every atom and committee member
 * straight from the scaled device descriptor G_dev_scaled and reads the
 * per-atom energies and descriptor gradients back to the host. The single
 * kernel below performs both the forward pass (the network energy of each atom)
 * and the backward pass (the gradient of that energy with respect to each input
 * descriptor), mirroring the host routines that evaluate the same feed-forward
 * network and its activation derivatives.
 */

/*******************************************************************************
 * \brief Atom-batched committee-network kernel: one thread block per
 *        (atom, committee) pair does a full forward pass followed by a backward
 *        pass. blockIdx.x selects the 0-based atom, blockIdx.y the 0-based
 *        committee member, and the block's threads split each layer's nodes by
 *        threadIdx.x stride.
 *
 *        Input and output layouts:
 *          - G_batch[i + n_input_max * atom]: scaled descriptor i of atom.
 *          - atomic_energy[c + n_committee * atom]: network energy of the atom
 *            for committee member c. The kernel adds atom_energies[ele] to the
 *            network output itself, so the result is the total atomic energy.
 *          - dE_dG[m + n_input_max * c + n_input_max * n_committee * atom]:
 *            gradient of that energy with respect to input descriptor m.
 *
 *        The weights are packed as
 *          W[i + max_n_nodes * j
 *            + max_n_nodes * max_n_nodes * c
 *            + max_n_nodes * max_n_nodes * n_committee * (L - 2)
 *            + max_n_nodes * max_n_nodes * n_committee * (n_layer - 1) * ele]
 *        and the biases as
 *          b[j + max_n_nodes * c
 *            + max_n_nodes * n_committee * (L - 2)
 *            + max_n_nodes * n_committee * (n_layer - 1) * ele],
 *        with n_nodes[(L - 1) + n_layer * ele] the node count of layer L for
 *        element ele. Activation codes match the enum in
 *        nnp_environment_types.F; normnodes divides each node's pre-activation
 *        (and its derivative) by the previous layer's node count.
 *
 *        Shared memory holds the per-layer pre- and post-activation node values
 *        and two alternating buffers for the chained descriptor gradient. The
 *        layout needs 2*max_n_nodes*n_layer + 2*n_input_max*max_n_nodes
 *        doubles per block; the launcher opts into the device's per-block
 *        shared-memory maximum when that exceeds the 48 KB default (e.g.
 *        80-input nets need ~105 KB, fine on an A100) and rejects only what
 *        the device cannot hold. The final layer is the scalar energy
 *        output (node 0).
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_nn_kernel_atomwise(
    const int n_atoms, const int n_committee, const int n_layer,
    const int max_n_nodes, const int n_input_max, const int normnodes,
    const int n_ele, const int *__restrict__ n_nodes,
    const int *__restrict__ act_fnct, const int *__restrict__ ele_ind,
    const double *__restrict__ atom_energies,
    const double *__restrict__ weights, const double *__restrict__ biases,
    const double *__restrict__ G_batch, double *__restrict__ atomic_energy,
    double *__restrict__ dE_dG) {

  const int atom = blockIdx.x;
  const int c = blockIdx.y;
  const int tid = threadIdx.x;
  const int nth = blockDim.x;
  if (atom >= n_atoms || c >= n_committee)
    return;

  const int ele = ele_ind[atom];
  const int n_input = n_nodes[0 + n_layer * ele];

  /* Shared memory layout. */
  extern __shared__ double smem[];
  double *node_pre = smem;
  double *node_post = node_pre + (size_t)max_n_nodes * n_layer;
  double *td_curr = node_post + (size_t)max_n_nodes * n_layer;
  double *td_prev = td_curr + (size_t)n_input_max * max_n_nodes;

  /* Init layer 1 (input descriptor) into node_post. 64-bit atom offset,
   * matching nnp_gpu_pack.h: G_batch scales with the rank's atom count. */
  for (int j = tid; j < n_input; j += nth) {
    node_post[j /* + 0*max_n_nodes */] =
        G_batch[j + (long long)n_input_max * atom];
  }
  __syncthreads();

  /* ---------- Forward ---------- */
  for (int L = 2; L <= n_layer; L++) {
    const int n_in = n_nodes[(L - 2) + n_layer * ele];
    const int n_out = n_nodes[(L - 1) + n_layer * ele];
    const size_t W_off =
        (size_t)max_n_nodes * max_n_nodes * c +
        (size_t)max_n_nodes * max_n_nodes * n_committee * (L - 2) +
        (size_t)max_n_nodes * max_n_nodes * n_committee * (n_layer - 1) * ele;
    const size_t b_off =
        (size_t)max_n_nodes * c + (size_t)max_n_nodes * n_committee * (L - 2) +
        (size_t)max_n_nodes * n_committee * (n_layer - 1) * ele;
    const double *W = weights + W_off;
    const double *b = biases + b_off;

    /* Each thread computes one output node j (stride nth). */
    for (int j = tid; j < n_out; j += nth) {
      double y = 0.0;
      for (int i = 0; i < n_in; i++) {
        y += W[i + max_n_nodes * j] * node_post[i + max_n_nodes * (L - 2)];
      }
      y += b[j];
      if (normnodes)
        y /= (double)n_in;

      /* Save pre-activation. */
      node_pre[j + max_n_nodes * (L - 1)] = y;

      /* Activation enum per nnp_environment_types.F. */
      const int act = act_fnct[L - 2];
      double y_post;
      switch (act) {
      case 1: /* tanh */
        y_post = tanh(y);
        break;
      case 2: /* gaus  f(x) = exp(-0.5 x^2) */
        y_post = exp(-0.5 * y * y);
        break;
      case 3: /* lin */
        y_post = y;
        break;
      case 4: /* cos */
        y_post = cos(y);
        break;
      case 5: /* sig  f(x) = 1 / (1 + exp(-x)) */
        y_post = 1.0 / (1.0 + exp(-y));
        break;
      case 6: /* invsig  f(x) = 1 - sig(x) */
        y_post = 1.0 - 1.0 / (1.0 + exp(-y));
        break;
      case 7: /* exp  f(x) = exp(-x) */
        y_post = exp(-y);
        break;
      case 8: /* softplus  f(x) = log(exp(x) + 1) */
        y_post = log(exp(y) + 1.0);
        break;
      case 9: /* quad  f(x) = x^2 */
        y_post = y * y;
        break;
      default:
        y_post = y;
        break;
      }
      node_post[j + max_n_nodes * (L - 1)] = y_post;
    }
    __syncthreads();
  }

  /* atomic_energy = node_post[layer n_layer, node 0] + atom_energies[ele]. */
  if (tid == 0) {
    atomic_energy[c + (long long)n_committee * atom] =
        node_post[0 + max_n_nodes * (n_layer - 1)] + atom_energies[ele];
  }

  /* ---------- Backward ---------- */
  /* Step 1: convert node_pre/node_post into activation derivative,
   * stored in node_pre (overwrite, pre-activation not needed
   * after this). Mirrors nnp_gradients_apply_actderiv. */
  for (int L = 2; L <= n_layer; L++) {
    const int n_out = n_nodes[(L - 1) + n_layer * ele];
    const int act = act_fnct[L - 2];
    for (int j = tid; j < n_out; j += nth) {
      const double pre = node_pre[j + max_n_nodes * (L - 1)];
      const double post = node_post[j + max_n_nodes * (L - 1)];
      double deriv;
      switch (act) {
      case 1: /* tanh: 1 - post^2 */
        deriv = 1.0 - post * post;
        break;
      case 2: /* gaus: -post * pre */
        deriv = -post * pre;
        break;
      case 3: /* lin: 1 */
        deriv = 1.0;
        break;
      case 4: /* cos: -sin(pre) */
        deriv = -sin(pre);
        break;
      case 5: /* sig: exp(-pre) / (1+exp(-pre))^2 */
      {
        const double e = exp(-pre);
        const double d = 1.0 + e;
        deriv = e / (d * d);
        break;
      }
      case 6: /* invsig: -sig'(pre) */
      {
        const double e = exp(-pre);
        const double d = 1.0 + e;
        deriv = -e / (d * d);
        break;
      }
      case 7: /* exp: -post */
        deriv = -post;
        break;
      case 8: /* softplus: (exp(post) + 1) / exp(post), uses POST per CPU code
               */
        deriv = (exp(post) + 1.0) / exp(post);
        break;
      case 9: /* quad: 2 * pre */
        deriv = 2.0 * pre;
        break;
      default:
        deriv = 1.0;
        break;
      }
      if (normnodes)
        deriv /= (double)n_nodes[(L - 2) + n_layer * ele];
      node_pre[j + max_n_nodes * (L - 1)] = deriv;
    }
  }
  __syncthreads();

  /* Step 2: init layer-2 tmp_der.
   *   tmp_der[2][i, j] = node_grad[2, j] * weights[2][i, j, c]
   * tmp_der_curr is used as the "current" buffer; it alternates
   * with tmp_der_prev each iteration of the chain. */
  {
    const int n_out_2 = n_nodes[1 + n_layer * ele];
    const size_t W2_off =
        (size_t)max_n_nodes * max_n_nodes * c +
        /* (L-2)=0 for L=2 */
        (size_t)max_n_nodes * max_n_nodes * n_committee * (n_layer - 1) * ele;
    const double *W2 = weights + W2_off;
    /* Each thread handles a row m of tmp_der. */
    for (int m = tid; m < n_input; m += nth) {
      for (int j = 0; j < n_out_2; j++) {
        td_curr[m + n_input_max * j] =
            node_pre[j + max_n_nodes * 1] /* deriv at layer 2 */
            * W2[m + max_n_nodes * j];
      }
    }
  }
  __syncthreads();

  /* Step 3: chain through layers 3..n_layer.
   *   tmp_der_new[m, j] = sum_kp tmp_der_old[m, kp] * weights[L][kp, j, c]
   *   then *= deriv[L, j] (Hadamard). */
  for (int L = 3; L <= n_layer; L++) {
    const int n_in_L = n_nodes[(L - 2) + n_layer * ele];
    const int n_out_L = n_nodes[(L - 1) + n_layer * ele];
    const size_t WL_off =
        (size_t)max_n_nodes * max_n_nodes * c +
        (size_t)max_n_nodes * max_n_nodes * n_committee * (L - 2) +
        (size_t)max_n_nodes * max_n_nodes * n_committee * (n_layer - 1) * ele;
    const double *WL = weights + WL_off;

    /* Swap roles: read td_curr (= old), write td_prev (= new). */
    for (int m = tid; m < n_input; m += nth) {
      for (int j = 0; j < n_out_L; j++) {
        double acc = 0.0;
        for (int kp = 0; kp < n_in_L; kp++) {
          acc += td_curr[m + n_input_max * kp] * WL[kp + max_n_nodes * j];
        }
        td_prev[m + n_input_max * j] =
            acc * node_pre[j + max_n_nodes * (L - 1)];
      }
    }
    __syncthreads();

    /* Swap pointers: td_curr <- td_prev (logical swap). */
    double *tmp_swap = td_curr;
    td_curr = td_prev;
    td_prev = tmp_swap;
  }

  /* Step 4: write dE_dG = tmp_der at last layer, j=0 column.
   *   For final layer with n_out = 1 (single energy output), tmp_der
   *   has shape (n_input, 1). dE_dG[m] = td_curr[m, 0]. */
  for (int m = tid; m < n_input; m += nth) {
    dE_dG[m + n_input_max * c + (long long)n_input_max * n_committee * atom] =
        td_curr[m + n_input_max * 0];
  }
}

/*******************************************************************************
 * \brief Once-at-init upload of the committee-network weights, biases,
 *        topology, activation functions and per-element atomic energies into
 *        the persistent state buffers (nn_n_nodes_dev, nn_act_fnct_dev,
 *        nn_atom_energies_dev, nn_W_packed_dev, nn_b_packed_dev), so the
 *        per-step launcher skips these host-to-device copies. Idempotent on the
 *        network shape: a second call with the same (n_committee, n_layer,
 *        max_n_nodes, n_ele) returns immediately, and a call with a different
 *        shape is an error. The packed layouts match the kernel's indexing
 *        (see nnp_nn_kernel_atomwise). The stream is synchronized before return
 *        so the caller's host staging arrays can go out of scope. Returns 0 on
 *        success, -1 on failure.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int
nnp_gpu_upload_nn_packed_c(int n_committee, int n_layer, int max_n_nodes,
                           int n_ele, const int *h_n_nodes,
                           const int *h_act_fnct, const double *h_atom_energies,
                           const double *h_W_packed, const double *h_b_packed) {
  if (n_committee <= 0 || n_layer < 2 || max_n_nodes <= 0 || n_ele <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_nn_packed: bad sizes "
            "n_committee=%d n_layer=%d max_n_nodes=%d n_ele=%d\n",
            n_committee, n_layer, max_n_nodes, n_ele);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  if (state->nn_packed_uploaded) {
    if (state->nn_packed_n_committee == n_committee &&
        state->nn_packed_n_layer == n_layer &&
        state->nn_packed_max_n_nodes == max_n_nodes &&
        state->nn_packed_n_ele == n_ele) {
      return 0;
    }
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_nn_packed: shape mismatch on second "
            "upload (n_committee %d->%d, n_layer %d->%d, "
            "max_n_nodes %d->%d, n_ele %d->%d)\n",
            state->nn_packed_n_committee, n_committee, state->nn_packed_n_layer,
            n_layer, state->nn_packed_max_n_nodes, max_n_nodes,
            state->nn_packed_n_ele, n_ele);
    return -1;
  }

  const size_t n_nodes_cnt = (size_t)n_layer * (size_t)n_ele;
  const size_t act_cnt = (size_t)(n_layer - 1);
  const size_t atom_e_cnt = (size_t)n_ele;
  const size_t W_cnt = (size_t)max_n_nodes * (size_t)max_n_nodes *
                       (size_t)n_committee * (size_t)(n_layer - 1) *
                       (size_t)n_ele;
  const size_t b_cnt = (size_t)max_n_nodes * (size_t)n_committee *
                       (size_t)(n_layer - 1) * (size_t)n_ele;

  if (nnp_gpu_pbuf_ensure((void **)&state->nn_n_nodes_dev,
                          &state->nn_n_nodes_dev_cap,
                          sizeof(int) * n_nodes_cnt))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->nn_act_fnct_dev,
                          &state->nn_act_fnct_dev_cap, sizeof(int) * act_cnt))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->nn_atom_energies_dev,
                          &state->nn_atom_energies_dev_cap,
                          sizeof(double) * atom_e_cnt))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->nn_W_packed_dev,
                          &state->nn_W_packed_dev_cap, sizeof(double) * W_cnt))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->nn_b_packed_dev,
                          &state->nn_b_packed_dev_cap, sizeof(double) * b_cnt))
    return -1;

  offloadMemcpyAsyncHtoD(state->nn_n_nodes_dev, h_n_nodes,
                         sizeof(int) * n_nodes_cnt, state->stream);
  offloadMemcpyAsyncHtoD(state->nn_act_fnct_dev, h_act_fnct,
                         sizeof(int) * act_cnt, state->stream);
  offloadMemcpyAsyncHtoD(state->nn_atom_energies_dev, h_atom_energies,
                         sizeof(double) * atom_e_cnt, state->stream);
  offloadMemcpyAsyncHtoD(state->nn_W_packed_dev, h_W_packed,
                         sizeof(double) * W_cnt, state->stream);
  offloadMemcpyAsyncHtoD(state->nn_b_packed_dev, h_b_packed,
                         sizeof(double) * b_cnt, state->stream);
  offloadStreamSynchronize(state->stream);

  state->nn_packed_uploaded = 1;
  state->nn_packed_n_committee = n_committee;
  state->nn_packed_n_layer = n_layer;
  state->nn_packed_max_n_nodes = max_n_nodes;
  state->nn_packed_n_ele = n_ele;
  return 0;
}

/*******************************************************************************
 * \brief Per-step committee-network launcher: evaluate every atom and committee
 *        member straight from the scaled device descriptor state->G_dev_scaled,
 *        then read the per-atom energies and descriptor gradients back to the
 *        host. The network weights, biases, topology, activation functions and
 *        atomic energies are read from the once-at-init device cache
 *        (state->nn_*_dev, populated by nnp_gpu_upload_nn_packed_c), so the
 *only per-step input that varies is the rank-local element index of each atom
 *(h_ele_ind), which is uploaded once and shared with the scaling launcher. The
 *grid is one block per (atom, committee) pair; the kernel writes
 *h_atomic_energy[c + n_committee * atom] and, when the caller passes a buffer
 *for it, h_dE_dG[m + n_input_max * c + n_input_max * n_committee * atom].
 *        h_dE_dG may be NULL, which leaves the gradient on the device for the
 *        force scatter and skips a read-back of several megabytes per step.
 *        The stream is synchronized after the read-backs so the caller can read
 *        the host arrays it asked for on return. Networks larger than the 48 KB
 *        default shared-memory window opt into the device per-block maximum;
 *        only sizes beyond that limit are rejected. Returns 0 on success, -1 on
 *        failure.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int
nnp_gpu_nn_atomwise_with_G_dev_c(const int n_atoms, const int n_committee,
                                 const int n_layer, const int max_n_nodes,
                                 const int n_input_max, const int normnodes,
                                 const int n_ele, const int *h_ele_ind,
                                 double *h_atomic_energy, double *h_dE_dG) {

  if (n_atoms <= 0 || n_committee <= 0 || n_layer < 2 || max_n_nodes <= 0 ||
      n_input_max <= 0 || n_ele <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_nn_atomwise_with_G_dev: bad sizes\n");
    return -1;
  }
  /* Shared-memory footprint of nnp_nn_kernel_atomwise for these sizes:
   * per-layer pre/post node values plus the two chained-gradient tiles. */
  const size_t smem_bytes =
      sizeof(double) * (2ULL * (size_t)max_n_nodes * (size_t)n_layer +
                        2ULL * (size_t)n_input_max * (size_t)max_n_nodes);

  /* Networks beyond the 48 KB default dynamic-shared-memory window (the
   * old 64x64 cap) opt into the device's per-block maximum instead of
   * being rejected: an 80-input net needs ~105 KB, which fits the A100's
   * 164 KB opt-in limit (Ada tops out near 99 KB, so such nets are still
   * rejected there with the size message below).
   * The granted size is a property of the network, not of the device, so what
   * is remembered is the number of bytes already opted into rather than a
   * boolean. A process that runs NNP environments in sequence can raise the
   * requirement, and launching against a smaller earlier grant would fail; a
   * request at or below the standing grant needs no new opt-in. */
  static size_t smem_optin_bytes = 0;
  if (smem_bytes > 48ULL * 1024ULL && smem_bytes > smem_optin_bytes) {
#if defined(__OFFLOAD_CUDA)
    int dev = 0, smem_max = 0;
    OFFLOAD_CHECK(cudaGetDevice(&dev));
    OFFLOAD_CHECK(cudaDeviceGetAttribute(
        &smem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev));
    if ((size_t)smem_max < smem_bytes ||
        cudaFuncSetAttribute(nnp_nn_kernel_atomwise,
                             cudaFuncAttributeMaxDynamicSharedMemorySize,
                             (int)smem_bytes) != cudaSuccess) {
      fprintf(stderr,
              "ERROR: nnp_gpu_nn_atomwise_with_G_dev: max_n_nodes=%d "
              "n_input_max=%d needs %zu bytes of shared memory per block; "
              "device allows %d.\n",
              max_n_nodes, n_input_max, smem_bytes, smem_max);
      return -1;
    }
#elif defined(__OFFLOAD_HIP)
    int dev = 0, smem_max = 0;
    OFFLOAD_CHECK(hipGetDevice(&dev));
    OFFLOAD_CHECK(hipDeviceGetAttribute(
        &smem_max, hipDeviceAttributeMaxSharedMemoryPerBlock, dev));
    if ((size_t)smem_max < smem_bytes) {
      fprintf(stderr,
              "ERROR: nnp_gpu_nn_atomwise_with_G_dev: max_n_nodes=%d "
              "n_input_max=%d needs %zu bytes of shared memory per block; "
              "device allows %d.\n",
              max_n_nodes, n_input_max, smem_bytes, smem_max);
      return -1;
    }
#endif
    smem_optin_bytes = smem_bytes;
  }

  const size_t E_cnt = (size_t)n_committee * (size_t)n_atoms;
  const size_t dEdG_cnt =
      (size_t)n_input_max * (size_t)n_committee * (size_t)n_atoms;

  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;
  if (state->G_dev_scaled == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_nn_atomwise_with_G_dev: G_dev_scaled "
                    "not initialised (call nnp_gpu_begin_descriptor_step + "
                    "nnp_gpu_scale_G_dev_c first)\n");
    return -1;
  }
  if (!state->nn_packed_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_nn_atomwise_with_G_dev: NN weights not "
                    "uploaded (call nnp_gpu_upload_nn_packed_c at NNP init)\n");
    return -1;
  }
  if (state->nn_packed_n_committee != n_committee ||
      state->nn_packed_n_layer != n_layer ||
      state->nn_packed_max_n_nodes != max_n_nodes ||
      state->nn_packed_n_ele != n_ele) {
    fprintf(stderr,
            "ERROR: nnp_gpu_nn_atomwise_with_G_dev: cached shape "
            "(n_committee=%d n_layer=%d max_n_nodes=%d n_ele=%d) "
            "differs from kernel call (%d/%d/%d/%d)\n",
            state->nn_packed_n_committee, state->nn_packed_n_layer,
            state->nn_packed_max_n_nodes, state->nn_packed_n_ele, n_committee,
            n_layer, max_n_nodes, n_ele);
    return -1;
  }

  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atomE_out,
                          &state->pbuf_atomE_out_cap, sizeof(double) * E_cnt))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_dEdG_out,
                          &state->pbuf_dEdG_out_cap, sizeof(double) * dEdG_cnt))
    return -1;
  /* The kernel writes only the rows below each element's own input width, so
   * for a model whose elements have different descriptor counts the rows in
   * between stay untouched. The force scatter never reads them, so the zeroing
   * is worth its bandwidth only when the caller is going to copy the whole
   * buffer back and would otherwise receive uninitialised device memory. */
  if (h_dE_dG != NULL)
    offloadMemsetAsync(state->pbuf_dEdG_out, 0, sizeof(double) * dEdG_cnt,
                       state->stream);

  /* Ensure the rank-local element index is on the device. Shared with the
   * scaling launcher, which applies the same guard, so only one host-to-device
   * copy happens per step across both. */
  if (nnp_gpu_step_ele_ind_ensure(state, h_ele_ind, n_atoms))
    return -1;

  /* Every kernel input except step_ele_ind_dev lives in state->nn_*_dev,
   * uploaded once by nnp_gpu_upload_nn_packed_c. */
  {
    dim3 grid(n_atoms, n_committee, 1);
    dim3 blk(32, 1, 1);
    nnp_nn_kernel_atomwise<<<grid, blk, smem_bytes, state->stream>>>(
        n_atoms, n_committee, n_layer, max_n_nodes, n_input_max, normnodes,
        n_ele, state->nn_n_nodes_dev, state->nn_act_fnct_dev,
        state->step_ele_ind_dev, state->nn_atom_energies_dev,
        state->nn_W_packed_dev, state->nn_b_packed_dev, state->G_dev_scaled,
        state->pbuf_atomE_out, state->pbuf_dEdG_out);
  }
  OFFLOAD_CHECK(offloadPeekAtLastError());

  offloadMemcpyAsyncDtoH(h_atomic_energy, state->pbuf_atomE_out,
                         sizeof(double) * E_cnt, state->stream);
  /* The descriptor gradient is optional: the force scatter reads the device
   * copy, so the whole-step driver passes NULL and skips a read-back that runs
   * to megabytes per step on the one synchronisation point of the step. */
  if (h_dE_dG != NULL)
    offloadMemcpyAsyncDtoH(h_dE_dG, state->pbuf_dEdG_out,
                           sizeof(double) * dEdG_cnt, state->stream);
  /* Sync after the D2Hs: the caller reads h_atomic_energy (and h_dE_dG when it
   * asked for it) on return, so the read-after-return contract requires the
   * copies to complete here. */
  offloadStreamSynchronize(state->stream);
  return 0;
}

#endif /* defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP) */
