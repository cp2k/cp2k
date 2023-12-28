# Nequip and Allegro

This
[Colab tutorial](https://colab.research.google.com/github/gabriele16/cp2k/blob/nequip-cp2k-colab/colab/allegro-cp2k-tutorial.ipynb)
illustrates how to train an equivariant neural network interatomic potential for bulk water using
the Allegro framework. You will learn how to train a model, deploy it in production, and run
molecular dynamics simulations in CP2K. The training and inference will be carried out on the GPU
provided by the Colab environment.

Allegro is designed for constructing highly accurate and scalable interatomic potentials for
molecular dynamics simulations. The methodology is described in detail in this paper
([](#Musaelian2023)). An open-source package that implements
[Allegro](https://github.com/mir-group/allegro), built on the
[Nequip framework](https://github.com/mir-group/nequip) was developed by the Allegro and NequIP
authors, A. Musaelian, S. Batzner, A. Johansson, L. Sun, C. J. Owen, M. Kornbluth, B. Kozinsky.

## Input Section

Inference in CP2K is performed through the
[NEQUIP](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.NEQUIP) and
[ALLEGRO](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.ALLEGRO) sections. As an example, the
relevant section for Allegro (or similarly for NequIP) is:

```none
&ALLEGRO
  ATOMS Si
  PARM_FILE_NAME Allegro/si-deployed.pth
  UNIT_COORDS angstrom
  UNIT_ENERGY eV
  UNIT_FORCES eV*angstrom^-1
&END ALLEGRO
```

where the `si-deployed.pth` refers to the PyTorch model that was deployed using the Allegro
framework, and the `UNIT` tags refer to the units of the coordinates, energy and forces of the model
itself. An example for the full input file can be found in the
[Colab tutorial](https://colab.research.google.com/github/gabriele16/cp2k/blob/nequip-cp2k-colab/colab/allegro-cp2k-tutorial.ipynb)
and on the regtests, see
[Allegro_si_MD.inp](https://github.com/cp2k/cp2k/blob/master/tests/Fist/regtest-allegro/Allegro_si_MD.inp)

### Input details

The tag [ATOMS](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.ALLEGRO.ATOMS) expects a list of
elements/kinds in a way and order that is consistent with the YAML file of NequIP and Allegro. If
this is not done unphysical results will be obtained. Additionally, the atomic coordinates in the
[COORD](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COORD) or [TOPOLOGY](#CP2K_INPUT.FORCE_EVAL.SUBSYS.TOPOLOGY)
section have to be provided in a way that is consistent with the YAML file. If this is not done
unphysical results will be obtained. Spotting such issues is quite straightforward as the energy is
significantly wrong. For example, by inverting the order of one of the elements in the test
`regtest-nequip/NequIP_water.inp`, the error with respect to the reference value is of the order of
1 eV. Additionally, running MD leads rapidly to highly unstable simulations.

## Compiling CP2K with LibTorch

Running with NequIP or Allegro requires compiling CP2K with the libtorch library. For the CP2K
binaries running on CPUs installing the toolchain using the flag `--with-libtorch` is enough. To
benefit from (often significant) GPU acceleration, the precompiled Libtorch library for CUDA can be
obtained at [https://pytorch.org](https://pytorch.org), for example for CUDA 11.8:

```shell
wget https://download.pytorch.org/libtorch/cu118/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcu118.zip
```

After extracting the libtorch CUDA binaries, the toolchain script `./install_cp2k_toolchain.sh` can
be run providing the appropriate path with the flag `--with-libtorch=<path-to-libtorch-cuda>`.

## Further Resources

For additional references on NequIP, Allegro and equivariant neural networks (e3nn) see:

- Allegro paper [](#Musaelian2023) and code <https://github.com/mir-group/allegro>
- NequIP paper [](#Batzner2022) and code <https://github.com/mir-group/nequip>
- A Tutorial on LAMMPS by the NequIP/Allegro authors is found at the Colab notebook
  [here](https://colab.research.google.com/drive/1yq2UwnET4loJYg_Fptt9kpklVaZvoHnq)
- For an introduction to e3nn see [e3nn.org](https://e3nn.org) and
  [doi:10.5281/zenodo.7430260](https://dx.doi.org/10.5281/zenodo.7430260)
