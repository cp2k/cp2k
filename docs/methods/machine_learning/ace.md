# Atomic Cluster Expansion (ACE)

The atomic cluster expansion is a complete descriptor of the local atomic environment. By
introducing nonlinear functions of the atomic cluster expansion an interatomic potential is obtained
that is comparable in accuracy to state-of-the-art machine learning potentials.

## Input Section

TODO

Inference in CP2K is performed through the [ACE](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.ACE)
section. As an example, the relevant section for ACE is:

```none
&ACE
  ATOMS O H
  POT_FILE_NAME ./sample.yaml
&END ACE
```

where the `sample.yaml` refers to the ACE model that was deployed using ACE. An example for the full
input file can be found in the regtests, see
[H2O-64_ACE_MD.inp](https://github.com/cp2k/cp2k/blob/master/tests/Fist/regtest-ace/H2O-64_ACE_MD.inp)

### Input details

The tag [ATOMS](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.ACE.ATOMS) expects a list of
elements/kinds that are to be treated with ACE.
[POT_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.ACE.POT_FILE_NAME) expects the
filename of the particular ACE fit.

## Compiling CP2K with ACE

Running with ACE requires compiling CP2K with the ACE library. For the CP2K binaries, please install
the toolchain using the flag `--with-ace`, which would download ACE from ACE Github release and
compile. GPU support is enabled when CUDA environment exists.

## Further Resources

For additional references on ACE see:

- ACE paper [](#Drautz2019), [](#Lysogorskiy2021), [](#Bochkarev2024) and code
  <https://github.com/ICAMS/lammps-user-pace>
