# How to compile the CP2K code

## 1. Acquire the code

For users, the preferred method is to [download a release](https://github.com/cp2k/cp2k/releases/)
(use the versioned tarballs, `cp2k-X.Y.tar.bz2`). For developers, the preferred method is to
[download from Git](./README.md#downloading-cp2k-source-code).

For more details on downloading CP2K, see <https://www.cp2k.org/download>.

## 2. Install prerequisites

The easiest way to build CP2K with all its dependencies is as a
[Docker container](./tools/docker/README.md).

CP2K supports GPU acceleration via CUDA (for NVIDIA GPUs), HIP/ROCm (for AMD GPUs), and OpenCL (for
a range of devices). If you wish to build with GPU support, please ensure to review sections 2i
(CUDA), 2w (ROCm/HIP), and 2x (OpenCL) in this document for detailed instructions.

Alternatively, the [toolchain script](./tools/toolchain/install_cp2k_toolchain.sh) can also be run
directly.

For a complete introduction to the toolchain script, see the [README](./tools/toolchain/README.md).

The basic steps are:

- Read toolchain installation options:

```shell
cd tools/toolchain/
./install_cp2k_toolchain.sh --help
```

- Launch toolchain script (example option choice for NVIDIA/CUDA):

```shell
./install_cp2k_toolchain.sh --with-libxsmm=install --with-openblas=system \
     --with-fftw=system --with-reflapack=no  --enable-cuda
```

- For AMD/ROCm/HIP support, use:

```shell
./install_cp2k_toolchain.sh --with-libxsmm=install --with-openblas=system \
     --with-fftw=system --with-reflapack=no  --enable-hip
```

## 3. Compile

This section has been moved to the
[manual](https://manual.cp2k.org/trunk/getting-started/build-from-source.html).

## 4. If it doesn't work

If things fail, take a break... go back to step 2.

## 5. Regtesting

If compilation works fine, it is recommended to test the generated binary, to exclude errors in
libraries, or miscompilations, etc.

```shell
./tests/do_regtest.py ./build/bin/ ${VERSION}
```

should work if you can locally execute CP2K without the need for, e.g., batch submission.

In the other case, you might need to configure the underlying testing script as described more
systematically at <https://www.cp2k.org/dev:regtesting>

## 6. Talk to us

In any case please tell us your comments, praise, criticism, thanks, etc. see
<https://www.cp2k.org>.

## 7. Manual

A reference manual of CP2K can be found on the web: <https://manual.cp2k.org> or can be generated
using the cp2k binary, see <https://manual.cp2k.org/trunk/generate_manual_howto.html>

## 8. Happy computing

The CP2K team.
