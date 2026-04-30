# GAPW_ACCURATE_XCINT default-readiness coverage

`GAPW_ACCURATE_XCINT` remains opt-in in this change. The tests below document the
initial coverage used to move towards making accurate XC integration the GAPW
default in a future change.

| Coverage area | Baseline / existing coverage | Accurate-XC coverage added here |
| --- | --- | --- |
| GAPW energy and forces | `regtest-acc-1/h2o.inp` | `regtest-acc-1/h2o-fine-1.inp` |
| Fine XC grid | regular-grid GAPW tests in `regtest-acc-1` | `h2o-fine-1.inp`, `h2o_f01_fine.inp`, `ft3_fine.inp` |
| NLCC | NLCC potentials in non-accurate tests | `regtest-acc-1/h2o-fine-nlcc-1.inp` |
| mGGA / tau | GAPW_XC TPSS tests in `regtest-gapw_xc` | `h2o-fine-tpss-1.inp`, `h2o-fine-tpss-stress-1.inp` |
| Analytical stress | regular GAPW/GAPW_XC energy coverage | `h2o-fine-tpss-stress-1.inp`, `Be_GAPW_XC_stress.inp` |
| Local XC energy density | regular local energy output | `regtest-acc-1/h2o-fine-local-energy-1.inp` |
| TDDFPT forces | `regtest-acc-5/h2o_f01.inp` | `regtest-acc-5/h2o_f01_fine.inp` |
| ADMM-GAPW TDDFPT response | `regtest-acc-5/ft3.inp` | `regtest-acc-5/ft3_fine.inp` |
| KG embedding | existing `regtest-kg` GPW cases | `H2-libxc-gapw.inp`, `H2-libxc-gapw_xc.inp` |
| KG atomic potential | `regtest-kg/H2_KG-1.inp` | `H2_KG-1-gapw.inp`, `H2_KG-1-gapw_xc.inp` |

The new tests intentionally compare explicit reference energies, force-debug
quantities, or debug-force consistency checks rather than relying only on
successful execution. This keeps the current default unchanged while making
future reference updates for a default flip easier to audit.

Remaining known gaps before a default change:

- `KG_METHOD` with GAPW/GAPW_XC still excludes `TNADD_METHOD EMBED_RI`.
- GAPW/ADMM local energy and stress densities on regular grids contain only the
  regular-grid / soft contribution, not hard one-center terms.
- A future default flip should run a broader reference-update pass comparing
  old default (`GAPW_ACCURATE_XCINT F`) and explicit accurate integration
  (`GAPW_ACCURATE_XCINT T`) for representative GAPW, GAPW_XC, ADMM-GAPW, mGGA,
  NLCC, stress, and TDDFPT-force cases.
