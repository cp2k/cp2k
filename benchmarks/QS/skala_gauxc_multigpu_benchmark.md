# GauXC/SKALA GPU benchmark helper

The helper script `skala_gauxc_multigpu_benchmark.py` runs one or more CP2K inputs with selected MPI
rank counts and writes a compact `summary.tsv`. It is intended for quick GauXC/SKALA energy and
force timing comparisons on GPU nodes.

Example:

```bash
python3 benchmarks/QS/skala_gauxc_multigpu_benchmark.py \
  --cp2k /path/to/cp2k.psmp \
  --input tests/QS/regtest-gauxc/H2_SKALA_ENERGY.inp \
  --input tests/QS/regtest-gauxc/H2_NATIVE_SKALA_GPW_FORCE_DEBUG.inp \
  --label h2_energy \
  --label h2_force \
  --ranks 1 2 \
  --outdir skala-gauxc-multigpu-run \
  --launcher-arg --bind-to \
  --launcher-arg none \
  --fill-fraction 0.95 \
  --atom-chunk-size 3 \
  --device-memory-cap 12000000000 \
  --distributed-torch 0 \
  --gradient-mpi-runtime 1
```

The output directory contains one CP2K output file per run and a combined `summary.tsv` with wall
time, final total FORCE_EVAL energy, peak GPU memory, and peak GPU utilization.

Relevant switches:

- `--distributed-torch 0` sets `GAUXC_ONEDFT_DISTRIBUTED_TORCH=0`, which keeps OneDFT Torch
  inference on rank 0 when GauXC supports this fallback.
- `--gradient-mpi-runtime 1` sets the experimental `GAUXC_GRADIENT_MPI_RUNTIME=1` override. This
  asks CP2K to use the normal GauXC MPI runtime for OneDFT nuclear gradients instead of the
  conservative replicated single-rank gradient runtime.
- `--fill-fraction`, `--atom-chunk-size`, and `--device-memory-cap` forward the corresponding GauXC
  runtime tuning environment variables.

The helper does not pin individual MPI ranks to separate GPUs. Use the MPI launcher or job scheduler
for rank-to-GPU binding when the benchmark requires strict affinity.
