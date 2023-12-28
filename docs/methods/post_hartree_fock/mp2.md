# Møller–Plesset Perturbation Theory

The MP2 method computes

$$
E^{(2)} = - \sum_{ijab} \frac{(ia|jb)[2(ia|jb)-(ib|ja)]}{\epsilon_{a}+\epsilon_{b}-\epsilon_{i}-\epsilon_{j}}
$$

with the GPW method as described in [](#DelBen2013) and [](#DelBen2012).

MP2 is controlled via the [WF_CORRELATION](#CP2K_INPUT.ATOM.METHOD.XC.WF_CORRELATION) section.

```none
&WF_CORRELATION
  METHOD  RI_MP2_GPW
  &WFC_GPW
    CUTOFF      300
    REL_CUTOFF  50
    EPS_FILTER  1.0E-12
    EPS_GRID    1.0E-8
  &END
  MEMORY    1200
  NUMBER_PROC  1
&END
```
