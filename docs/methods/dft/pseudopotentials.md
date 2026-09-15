# Pseudopotentials

Most GPW calculations in CP2K use norm-conserving Goedecker-Teter-Hutter ({term}`GTH`)
pseudopotentials. A pseudopotential removes chemically inactive core electrons from the explicit
electronic problem and represents their effect on the valence electrons through an effective
potential. This reduces the number of electrons and avoids the very hard core density that would
otherwise require extremely fine grids.

Pseudopotential files are selected in
[POTENTIAL_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME), and the actual potential is
selected for each atomic [KIND](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND) with
[POTENTIAL](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.POTENTIAL):

```text
&FORCE_EVAL
  &DFT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
  &END DFT
  &SUBSYS
    &KIND O
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

Do not confuse the [POTENTIAL](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.POTENTIAL) *keyword* with the
identically-named [POTENTIAL](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.POTENTIAL_SECTION) *section* that
also defines a pseudopotential. Their distinction is that the *keyword* takes the type and name of
the potential, while the *section* takes the complete specification and data in an internal format.
Specifying only the *keyword* as is done above suffices for most practical usage. For those curious,
CP2K will resolve user-provided inputs and compose restart files with more verbose but otherwise
equivalent syntax as a new input (see [MOTION/PRINT/RESTART](#CP2K_INPUT.MOTION.PRINT.RESTART)), and
such a restart file gives a glimpse of the usage of both the *keyword* and the *section*:

```
     &KIND "O"
       POTENTIAL "GTH-PBE-q6"
       &POTENTIAL
         2 4
           2.4455430000000000E-001 2 -1.6667214800000000E+001  2.4873113199999999E+000
         2
           2.2095592000000000E-001 1  1.8337458110000000E+001
           2.1133246999999999E-001 0
         # Potential name: GTH-PBE-Q6 for element symbol: O
         # Potential read from the potential filename: GTH_POTENTIALS
       &END POTENTIAL
     &END KIND
     &KIND "H"
       POTENTIAL "GTH-PBE-q1"
       &POTENTIAL
         1
           2.0000000000000001E-001 2 -4.1789004399999996E+000  7.2446330999999997E-001
         0
         # Potential name: GTH-PBE-Q1 for element symbol: H
         # Potential read from the potential filename: GTH_POTENTIALS
       &END POTENTIAL
     &END KIND
```

The suffix `q6` in `GTH-PBE-q6`, for example, means that six valence electrons are treated
explicitly. The chosen basis set should match this valence configuration; for oxygen, the common
basis sets with the same suffix in the full name like `DZVP-MOLOPT-GTH-q6` can be used. If some
pseudopotentials without corresponding basis sets, or vice versa, are spotted in the built-in data
files, consult developers for help.

## Choosing a Pseudopotential

Use a pseudopotential generated for the exchange-correlation functional family used in the
calculation. For example, `GTH-PBE-q6` is a natural choice for PBE calculations with oxygen. Mixing
functional families can be acceptable for exploratory work in some cases, but it is not a systematic
route to high accuracy.

The CP2K data directory contains several pseudopotential libraries:

- `GTH_POTENTIALS` contains widely used GTH potentials for common GPW calculations.
- `POTENTIAL_UZH` contains the UZH protocol GTH potentials designed to be used with matching UZH
  basis sets.
- `NLCC_POTENTIALS` and `GTH_SOC_POTENTIALS` contain more specialized potentials.
- `ECP_POTENTIALS` contains effective core potentials for Gaussian integral based calculations.

For new GPW production inputs, prefer a matching UZH protocol pair from `POTENTIAL_UZH` and
`BASIS_MOLOPT_UZH` when it is available for the element and functional family. The older
`GTH_POTENTIALS` library remains important for reproducing established calculations and for cases
where a matching UZH setup is not available.

For all-electron calculations, use `POTENTIAL ALL` together with an all-electron basis set and the
GAPW method:

```text
&KIND O
  BASIS_SET SVP-MOLOPT-GGA-ae
  POTENTIAL ALL
&END KIND
```

## Consistency Checks

Useful checks when setting up a calculation are:

- The basis set and pseudopotential should be available in the files named in the `DFT` section.
- The pseudopotential valence charge should match the basis set suffix where such a suffix is used.
- The exchange-correlation functional should be consistent with the pseudopotential family.
- For heavy elements, decide whether a large-core, medium-core, small-core, or all-electron
  description is appropriate for the property of interest.

For a tested minimal GPW input using GTH pseudopotentials, see
[](../../getting-started/first-calculation).

## See Also

- <https://en.wikipedia.org/wiki/Pseudopotential>
- <https://cp2k.org/static/potentials/>
- <https://www.cp2k.org/tools:cp2k-basis>
- [](#Goedecker1996)
- [](#Hartwigsen1998)
- [](#Krack2005)
- [](#Iannuzzi2026)
