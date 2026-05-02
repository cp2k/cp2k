# GAPW_ACCURATE_XCINT default-readiness coverage

`GAPW_ACCURATE_XCINT` remains opt-in in this change. The tests below document the initial coverage
used to move towards making accurate XC integration the GAPW default in a future change.

| Coverage area              | Baseline / existing coverage                                                  | Accurate-XC coverage added here                                                                                      |
| -------------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| GAPW energy and forces     | `regtest-acc-1/h2o.inp`                                                       | `regtest-acc-1/h2o-fine-1.inp`, `h2o-fine-tpss-force-1.inp`                                                          |
| Fine XC grid               | regular-grid GAPW tests in `regtest-acc-1`                                    | `h2o-fine-1.inp`, `h2o_f01_fine.inp`, `ft3_fine.inp`                                                                 |
| nonlocal vdW               | rVV10 and vdW-DF/optB88 coverage in `regtest-dft-vdw-corr-*`                  | `regtest-acc-1/argon-rvv10-gapw-accurate.inp`, `argon-vdW-DF-optB88-gapw-accurate.inp`                               |
| NLCC                       | NLCC potentials in non-accurate tests                                         | `regtest-acc-1/h2o-fine-nlcc-1.inp`                                                                                  |
| mGGA / tau                 | GAPW_XC TPSS tests in `regtest-gapw_xc`                                       | `h2o-fine-tpss-1.inp`, `h2o-fine-tpss-force-1.inp`, `h2o-fine-tpss-stress-1.inp`, `h2o-fine-tpss-stress-debug-1.inp` |
| GAPW_XC                    | regular GAPW_XC tests in `regtest-gapw_xc`, including `Be_GAPW_XC_stress.inp` | `regtest-acc-1/Ar-2.inp`, `Ar-4.inp`, `h2o-gapw_xc-force-1.inp`, `h2o-gapw_xc-stress-debug-1.inp`                    |
| spin/open shell            | triplet GAPW tests in `regtest-acc-3`                                         | existing H2O/Ne triplets plus open-shell XAS_TDP and ADMM-GAPW O2 stress coverage                                    |
| k-points                   | GAPW/GAPW_XC k-point tests in `regtest-kp-1` and `regtest-kp-2`               | `regtest-acc-1/Ar-kpoints-gapw-accurate.inp`, `C-kpoints-gapw_xc-accurate.inp`                                       |
| Analytical/full stress     | regular GAPW/GAPW_XC energy coverage, including `Be_GAPW_XC_stress.inp`       | diagonal and full-matrix stress checks for GAPW, GAPW_XC, Fine-XC/mGGA, ADMM-GAPW/Fine-XC, and periodic/off-diagonal cases |
| Local XC energy density    | regular local energy output                                                   | `regtest-acc-1/h2o-fine-local-energy-1.inp`                                                                          |
| ADMM-GAPW stress FD check  | ADMM stress output in `regtest-admm-qps-2`                                    | diagonal stress checks plus `h2o-admm-gapw-pbe-fine-full-stress-debug-1.inp`                                         |
| ADMM-GAPW open shell       | ADMM-GAPW O2 stress output in `regtest-admm-qps-2`                            | `regtest-acc-2/O2-admmq-gapw-open-shell-accurate-stress.inp`                                                         |
| DC-DFT/Energy Correction   | `regtest-acc-2/HF-ec1.inp`, `HF-ec2.inp`, `HF-ec3.inp`                       | `HF-ec4.inp`, `HF-ec5.inp`, `HF-ec6.inp`, and `HF-ec7.inp` now also use explicit accurate integration                |
| TDDFPT forces              | `regtest-acc-5/h2o_f01.inp`                                                   | `regtest-acc-5/h2o_f01_fine.inp`                                                                                     |
| ADMM-GAPW TDDFPT response  | `regtest-acc-5/ft3.inp`                                                       | `regtest-acc-5/ft3_fine.inp`                                                                                         |
| XAS/RT response            | XAS_TDP and RTBSE coverage in `regtest-xastdp` and `regtest-rtbse`            | open-shell GAPW XAS_TDP and GAPW RTBSE smoke tests in `regtest-acc-3`                                                |
| KG embedding               | existing `regtest-kg` GPW cases                                               | GAPW/GAPW_XC energy and stress checks for libxc KG and RI embedding                                                  |
| KG atomic potential        | `regtest-kg/H2_KG-1.inp`                                                      | GAPW/GAPW_XC energy and stress checks for `TNADD_METHOD ATOMIC`                                                      |

The new tests intentionally compare explicit reference energies, force-debug quantities, or
debug-force/stress consistency checks rather than relying only on successful execution. This keeps
the current default unchanged while making future reference updates for a default flip easier to
audit.

Local one-rank smoke comparisons between the current default (`GAPW_ACCURATE_XCINT F`) and explicit
accurate integration (`GAPW_ACCURATE_XCINT T`) were run for representative cases:

| Case                     | Representative input                                    | F energy [Ha]      | T energy [Ha]      | T-F [Ha]          |
| ------------------------ | ------------------------------------------------------- | ------------------ | ------------------ | ----------------- |
| GAPW PBE                 | `regtest-acc-1/h2o.inp`                                 | -17.26123610393284 | -17.26044915752598 | 0.00078694640686 |
| GAPW_XC Pade             | `regtest-acc-1/Ar-2.inp`                                | -21.04944269549961 | -21.04944239056048 | 0.00000030493913 |
| mGGA TPSS Fine-XC        | `regtest-acc-1/h2o-fine-tpss-1.inp`                     | -17.26821914033707 | -17.26735923212099 | 0.00085990821608 |
| NLCC Fine-XC             | `regtest-acc-1/h2o-fine-nlcc-1.inp`                     | -18.07356942938059 | -18.07338831671611 | 0.00018111266448 |
| TPSS stress FD           | `regtest-acc-1/h2o-fine-tpss-stress-debug-1.inp`        | -17.17596758695775 | -17.17509714444003 | 0.00087044251771 |
| ADMM-GAPW PBE stress FD  | `regtest-acc-2/h2o-admm-gapw-pbe-stress-debug-1.inp`    | -21.03576759790435 | -21.11980618566643 | -0.08403858776208 |
| TDDFPT force Fine-XC     | `regtest-acc-5/h2o_f01_fine.inp`                        | -16.86015219162475 | -16.86018894515963 | -0.00003675353488 |
| KG RI GAPW               | `regtest-kg/H2_H2O-kglri-gapw.inp`                      | -18.31547980084748 | -18.31452827588301 | 0.00095152496447 |

Finite-difference checks with `STOP_ON_MISMATCH` are included for:

- TPSS/mGGA forces: `regtest-acc-1/HF-d2.inp`, `regtest-acc-1/h2o-fine-tpss-force-1.inp`.
- TPSS/mGGA diagonal stress: `regtest-acc-1/h2o-fine-tpss-stress-debug-1.inp`.
- Full analytical stress matrices with off-diagonal virial components: `regtest-acc-1/h2o-gapw-pbe-full-stress-debug-1.inp`,
  `h2o-fine-tpss-full-stress-debug-1.inp`, `h2o-gapw_xc-full-stress-debug-1.inp`,
  `regtest-acc-2/h2o-admm-gapw-pbe-fine-full-stress-debug-1.inp`.
- GAPW_XC forces and diagonal stress: `regtest-acc-1/h2o-gapw_xc-force-1.inp`,
  `regtest-acc-1/h2o-gapw_xc-stress-debug-1.inp`.
- ADMM-GAPW forces: `regtest-acc-1/HF-d5.inp`, `regtest-acc-5/ft3_fine.inp`.
- ADMM-GAPW diagonal stress: `regtest-acc-2/h2o-admm-gapw-stress-debug-1.inp`,
  `regtest-acc-2/h2o-admm-gapw-pbe-stress-debug-1.inp`.
- Fine-XC forces and diagonal stress: `regtest-acc-1/h2o-fine-tpss-force-1.inp`,
  `regtest-acc-1/h2o-fine-tpss-stress-debug-1.inp`, `regtest-acc-5/ft3_fine.inp`.
- KG GAPW/GAPW_XC stress smoke checks: `regtest-kg/H2-libxc-gapw-stress.inp`,
  `H2-libxc-gapw_xc-stress.inp`, `H2_KG-1-gapw-stress.inp`, `H2_KG-1-gapw_xc-stress.inp`,
  `H2_H2O-kglri-gapw-stress.inp`, `H2_H2O-kglri-gapw_xc-stress.inp`.
- DC-DFT/Energy Correction force-debug checks: `regtest-acc-2/HF-ec1.inp` through
  `regtest-acc-2/HF-ec7.inp`.

Remaining known gaps before a default change:

- A future default flip should still run the full CP2K regtest suite and update affected default
  references deliberately; the representative F/T checks above are the local default-readiness
  smoke matrix for this PR.
- SIC is not part of the green GAPW accurate-integration matrix: the implementation explicitly
  aborts regular SIC with GAPW (`sic and GAPW not yet compatible`), and GAPW_XC XC-potential paths
  assert `sic_none`.
- Low-spin ROKS with GAPW/GAPW_XC, and with ADMM, remains explicitly guarded as incompatible in
  `qs_ks_utils.F`.
- KG combinations that are guarded by the implementation remain out of scope for this coverage
  matrix: KG meta-kinetic energy functionals abort as not implemented, and GAPW/GAPW_XC with the
  Energy Correction/Harris response machinery still has independent GAPW guards.
- The added XAS_TDP and RTBSE tests are representative smoke coverage for GAPW accurate integration,
  not a full spectroscopy/response-method matrix over all kernels, SOC, GW/BSE, periodic, and
  open-shell variants.

Related print-key limitation that is independent of the `GAPW_ACCURATE_XCINT` default:

- `LOCAL_ENERGY_CUBE` and `LOCAL_STRESS_CUBE` are regular-grid print keys. For GAPW/GAPW_XC and
  ADMM-GAPW they keep the existing soft-grid semantics; atom-centered hard one-center terms are not
  projected onto the cube grid.
