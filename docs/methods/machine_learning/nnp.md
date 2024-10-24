# Neural Network Potentials

## Input Section

NNPs in CP2K are supported via the `NNP` Force_Eval section. As an example, a basic `NNP` section would be:

```
&NNP
  NNP_INPUT_FILE_NAME PATH/TO/input.nn
  SCALE_FILE_NAME PATH/TO/scaling.data
  &MODEL
    WEIGHTS PATH/TO/weights
  &END MODEL
&END NNP
```

Where `PATH/TO/` points to the relevant `NNP` input files.

In the meantime, the following links might be helpful:

- <https://www.cp2k.org/tools:aml>
- <https://doi.org/10.1063/5.0160326>
- [](#Behler2007)
- [](#Behler2011)
- [](#Schran2020)
- [](#Schran2020b)
